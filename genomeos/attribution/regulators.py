# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""M1b: which factor acts at each compiled element, so a compiled rule can name a gene as its source.

M1 (genomeos/attribution/runrules.py, data/results/runrules_census.json) found that 0 of 440,589
compiled rules names a gene as its source: every rule is element -> gene, so a run of the compiled
programme is one independent one-gene system per regulated gene and cannot show a network effect. A
regulator -> element -> gene rule needs the factor that binds the element. This module assigns, to
each compiled rule, the candidate factors with the kind of evidence behind each, and it never writes a
compiled program: `genomeos/attribution/compile.py` is the only writer of those, and what it must
change is stated in `COMPILE_SPEC`.

Evidence kinds, strongest first, and what each is not:

- `occupancy_same_cell_with_motif`: an ENCODE TF ChIP-seq peak cluster of the factor overlaps the
  element in the rule's own cell type, and the factor's own JASPAR profile matches inside the element.
  Measured occupancy, with a sequence match that makes direct binding plausible. Occupancy is not
  function: a bound factor need not change the gene.
- `occupancy_same_cell`: the same peak overlap, with no match of the factor's profile (or no profile).
  Measured occupancy, possibly indirect.
- `occupancy_other_cell`: the factor occupies the element only in cell types other than the rule's.
  Not a source in the rule's cell; recorded so it is not lost, never emitted as a rule there.
- `motif_only`: a JASPAR match and no occupancy anywhere. Predicted evidence: a motif is not occupancy.
- a rule with none of these keeps its source UNKNOWN.

Direction. The direction of a factor's effect on the gene needs its own evidence and is never defaulted
to activation. The only direction evidence used here is curated, factor-level and not element-level:
the factor's Gene Ontology molecular-function annotation as a DNA-binding transcription activator or
repressor. A factor -> gene rule gets a direction only when the factor has exactly one of the two and it
agrees with the element's own compiled rule (an activator at an element whose compiled rule activates
the gene; a repressor at one whose compiled rule inhibits it). Every other case is `unknown`, with the
reason.

The network-level test (registered in data/results/regulators_registration.json before its endpoint is
read): for a factor knocked down in a cell type, do the genes whose compiled rules in that cell pass
through elements the factor occupies there change more than the other genes of the same compiled
programme, more than genes near the factor's peaks (a baseline that ignores the compiled rules), and
more than another factor's assigned genes do (the matched control)? The functions below the
`# --- the registered test` line are pure: inputs are arguments, nothing is read from disk.
"""

from __future__ import annotations

import bisect
import gzip
import math
import random
import re
from collections import Counter, defaultdict
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

LANE = "lane-regulators"
CHROMS = tuple(f"chr{c}" for c in [*range(1, 23), "X", "Y"])
KNOWLEDGE = Path("data/knowledge/regulators")
CACHE = Path("data/cache/regulators")
CLUSTERS = KNOWLEDGE / "encRegTfbsClusteredWithCells.hg38.bed.gz"
KD_DIR = KNOWLEDGE / "knockdown"
GOA = Path("data/knowledge/goa_human.gaf.gz")

CLUSTER_SOURCE = {
    "name": "ENCODE 3 TF ChIP-seq peak clusters with cell types (UCSC track encRegTfbsClustered, hg38)",
    "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/encRegTfbsClustered/"
    "encRegTfbsClusteredWithCells.hg38.bed.gz",
    "track_page": "https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=encRegTfbsClustered",
    "bytes": 125_375_233,
    "sha256": "542eb15d4d90236c2d3437b285396bbaddee023ca3d9d0e7efb767fcd1ef7899",
    "last_modified": "2020-01-27",
    "assembly": "GRCh38",
    "licence": "ENCODE data, no restrictions on use (ENCODE data use policy); UCSC clustering of the "
    "ENCODE uniform peaks",
    "columns": "chrom, start, end, factor, score (0-1000), comma-separated ENCODE biosample names",
    "kind": "measured occupancy (ChIP-seq peaks); not function",
}

GO_ACTIVATOR = frozenset({"GO:0001228", "GO:0001216"})  # DNA-binding transcription activator activity
GO_REPRESSOR = frozenset({"GO:0001227", "GO:0001217"})  # DNA-binding transcription repressor activity
GO_SOURCE = {
    "file": str(GOA),
    "terms_activator": sorted(GO_ACTIVATOR),
    "terms_repressor": sorted(GO_REPRESSOR),
    "qualifiers_used": "enables and contributes_to; any NOT qualifier is skipped",
    "licence": "Gene Ontology annotations, CC BY 4.0",
    "kind": "curated, factor-level: says what the factor does at the genes it regulates in general, "
    "not at this element",
}

KINDS = (
    "occupancy_same_cell_with_motif",
    "occupancy_same_cell",
    "occupancy_other_cell",
    "motif_only",
)
SOURCE_KINDS = KINDS[:2]  # the two that give a rule a source in its own cell
CELL_ALIASES = {"H1-hESC": "H1"}  # ENCODE's older name for the biosample the compiled programs call H1


def context_label(encode_cell: str) -> str:
    """An ENCODE biosample name as the compiled programs spell a cell context (`HeLa-S3` -> `HeLa_S3`)."""
    name = CELL_ALIASES.get(encode_cell, encode_cell)
    return re.sub(r"[^A-Za-z0-9]", "_", name)


# --- compiled rules ------------------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class CompiledRule:
    program: str
    element: str
    gene: str
    action: str  # "activates" | "inhibits"
    cell: str  # the rule's `when: cell_type`, "" when it has none
    evidence: str  # the rule's evidence kind: "predicted" | "experimental"
    chrom: str
    start: int
    end: int


def compiled_rules(module: Any, program: str) -> list[CompiledRule]:
    """Every activates/inhibits rule of a parsed compiled program, with its element's locus."""
    out = []
    for r in module.rules:
        act = getattr(r.action, "value", str(r.action))
        if act not in ("activates", "inhibits"):
            continue
        e = module.entities.get(r.source)
        loc = getattr(e, "locus", None)
        if loc is None:
            raise ValueError(f"{program}: rule source {r.source!r} has no locus")
        out.append(
            CompiledRule(
                program,
                r.source,
                r.target,
                act,
                r.when.get("cell_type", "") if r.when else "",
                getattr(r.evidence.kind, "value", str(r.evidence.kind)),
                loc.chrom,
                int(loc.start),
                int(loc.end),
            )
        )
    return out


# --- measured occupancy: the ENCODE peak clusters ---------------------------------------------------------


@dataclass
class Clusters:
    """Peak clusters per chromosome, sorted by start, with a running maximum of ends for overlap."""

    factors: list[str]
    cellsets: list[frozenset[str]]  # compiled context labels
    starts: dict[str, list[int]] = field(default_factory=dict)
    ends: dict[str, list[int]] = field(default_factory=dict)
    factor_of: dict[str, list[int]] = field(default_factory=dict)
    cellset_of: dict[str, list[int]] = field(default_factory=dict)
    max_end: dict[str, list[int]] = field(default_factory=dict)

    def overlapping(self, chrom: str, start: int, end: int) -> list[int]:
        """Indices of clusters on `chrom` sharing at least one base with [start, end)."""
        st = self.starts.get(chrom)
        if not st:
            return []
        en, mx = self.ends[chrom], self.max_end[chrom]
        j = bisect.bisect_left(st, end) - 1
        out = []
        while j >= 0 and mx[j] > start:
            if en[j] > start:
                out.append(j)
            j -= 1
        return out

    def cells_of(self) -> Counter:
        """Clusters per compiled cell label (a cluster counts once per cell type it names)."""
        n: Counter = Counter()
        for chrom in self.cellset_of:
            for k, c in Counter(self.cellset_of[chrom]).items():
                for cell in self.cellsets[k]:
                    n[cell] += c
        return n

    def factors_in(self) -> dict[str, set[str]]:
        """Compiled cell label -> the factors with at least one cluster naming it."""
        out: dict[str, set[str]] = defaultdict(set)
        for chrom in self.cellset_of:
            for fi, ci in set(zip(self.factor_of[chrom], self.cellset_of[chrom], strict=True)):
                for cell in self.cellsets[ci]:
                    out[cell].add(self.factors[fi])
        return out


def parse_clusters(lines: Iterable[str], chroms: Iterable[str] = CHROMS) -> Clusters:
    keep = set(chroms)
    factors: dict[str, int] = {}
    cellsets: dict[str, int] = {}
    rows: dict[str, list[tuple[int, int, int, int]]] = defaultdict(list)
    for line in lines:
        f = line.rstrip("\n").split("\t")
        if len(f) < 6 or f[0] not in keep:
            continue
        fi = factors.setdefault(f[3], len(factors))
        ci = cellsets.setdefault(f[5], len(cellsets))
        rows[f[0]].append((int(f[1]), int(f[2]), fi, ci))
    sets = [frozenset(context_label(c) for c in k.split(",") if c) for k in cellsets]
    cl = Clusters(list(factors), sets)
    for chrom, rs in rows.items():
        rs.sort()
        cl.starts[chrom] = [r[0] for r in rs]
        cl.ends[chrom] = [r[1] for r in rs]
        cl.factor_of[chrom] = [r[2] for r in rs]
        cl.cellset_of[chrom] = [r[3] for r in rs]
        mx, run = [], -1
        for r in rs:
            run = max(run, r[1])
            mx.append(run)
        cl.max_end[chrom] = mx
    return cl


def load_clusters(path: Path = CLUSTERS, chroms: Iterable[str] = CHROMS) -> Clusters:
    with gzip.open(path, "rt") as fh:
        return parse_clusters(fh, chroms)


def element_occupancy(cl: Clusters, chrom: str, start: int, end: int) -> dict[str, frozenset[str]]:
    """Factor -> the compiled cell labels in which a cluster of that factor overlaps the element."""
    out: dict[str, set[str]] = defaultdict(set)
    for j in cl.overlapping(chrom, start, end):
        out[cl.factors[cl.factor_of[chrom][j]]] |= cl.cellsets[cl.cellset_of[chrom][j]]
    return {f: frozenset(c) for f, c in out.items()}


def factor_peaks(cl: Clusters, factor: str, cell: str) -> dict[str, list[int]]:
    """Midpoints of the factor's clusters that name `cell`, per chromosome, sorted."""
    fi = cl.factors.index(factor) if factor in cl.factors else -1
    out: dict[str, list[int]] = {}
    if fi < 0:
        return out
    for chrom, fs in cl.factor_of.items():
        cs, st, en = cl.cellset_of[chrom], cl.starts[chrom], cl.ends[chrom]
        mids = [(st[j] + en[j]) // 2 for j, f in enumerate(fs) if f == fi and cell in cl.cellsets[cs[j]]]
        if mids:
            out[chrom] = sorted(mids)
    return out


# --- predicted: JASPAR matches -----------------------------------------------------------------------


def motif_factor_names(jaspar_name: str) -> list[str]:
    """The HGNC-style factor names one JASPAR profile name stands for.

    `MAX(var.2)` -> [MAX]; a heterodimer `GATA1::TAL1` -> [GATA1, TAL1]: a match of a dimer's profile is
    credited to each partner, because either can be the one that occupies the site."""
    base = re.sub(r"\(VAR\.\d+\)", "", jaspar_name.upper()).strip()
    return [p.strip() for p in base.split("::") if p.strip()]


# --- direction: curated, factor-level ------------------------------------------------------------------


def parse_go_directions(lines: Iterable[str]) -> dict[str, dict[str, Any]]:
    """Symbol -> {activator, repressor, evidence codes} from a GAF 2.2 stream."""
    out: dict[str, dict[str, Any]] = {}
    for line in lines:
        if line.startswith("!"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 7 or "NOT" in f[3]:
            continue
        go = f[4]
        if go not in GO_ACTIVATOR and go not in GO_REPRESSOR:
            continue
        d = out.setdefault(f[2], {"activator": False, "repressor": False, "codes": set()})
        d["activator" if go in GO_ACTIVATOR else "repressor"] = True
        d["codes"].add(f[6])
    return out


def load_go_directions(path: Path = GOA) -> dict[str, dict[str, Any]]:
    with gzip.open(path, "rt") as fh:
        return parse_go_directions(fh)


def factor_class(go: dict[str, dict[str, Any]], factor: str) -> str:
    d = go.get(factor)
    if not d:
        return "no_curated_direction"
    if d["activator"] and d["repressor"]:
        return "activator_and_repressor"
    return "activator" if d["activator"] else "repressor"


def rule_direction(fclass: str, element_action: str) -> tuple[str, str]:
    """The direction of a factor -> gene rule through an element, and why; never a default."""
    if fclass == "activator" and element_action == "activates":
        return "activates", "curated activator at an element whose compiled rule activates the gene"
    if fclass == "repressor" and element_action == "inhibits":
        return "inhibits", "curated repressor at an element whose compiled rule inhibits the gene"
    if fclass in ("activator", "repressor"):
        return "unknown", f"curated {fclass} at an element whose compiled rule {element_action} the gene"
    return "unknown", fclass


# --- the assignment and the census ---------------------------------------------------------------------


def candidates(
    rule: CompiledRule,
    occupancy: dict[str, frozenset[str]],
    motif: frozenset[str],
    go: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    """The factors that occupy a rule's element, each with its evidence kind relative to the rule's cell.

    Only occupied factors are listed one by one; motif-only matches are counted by the census, not
    listed, because a match alone names hundreds of factors per element."""
    out = []
    for f in sorted(occupancy):
        cells = occupancy[f]
        own_cell = bool(rule.cell) and rule.cell in cells
        kind = (KINDS[0] if f in motif else KINDS[1]) if own_cell else KINDS[2]
        direction, basis = (
            rule_direction(factor_class(go, f), rule.action) if kind in SOURCE_KINDS else ("unknown", "")
        )
        out.append(
            {
                "factor": f,
                "kind": kind,
                "motif": f in motif,
                "occupied_in": len(cells),
                "direction": direction,
                "direction_basis": basis,
            }
        )
    return out


def rule_kind(cands: list[dict[str, Any]], motif_count: int) -> str:
    """The strongest evidence kind any candidate gives the rule relative to its own cell, or `unknown`."""
    kinds = {c["kind"] for c in cands}
    for k in KINDS:
        if k in kinds:
            return k
    return "motif_only" if motif_count > 0 else "unknown"


def empty_tally() -> dict[str, Any]:
    return {
        "rules": 0,
        "rules_gaining_a_source": 0,
        "by_kind": {k: 0 for k in (*KINDS, "unknown")},
        "candidate_pairs_in_own_cell": 0,
        "candidate_pairs_with_motif": 0,
        "direction": Counter(),
        "rules_with_a_directed_candidate": 0,
        "factors": Counter(),
        "sources_per_rule": [],
        "motif_factors_per_rule": [],
    }


def tally(t: dict[str, Any], rule: CompiledRule, cands: list[dict[str, Any]], motif_count: int) -> None:
    own = [c for c in cands if c["kind"] in SOURCE_KINDS]
    t["rules"] += 1
    t["by_kind"][rule_kind(cands, motif_count)] += 1
    t["motif_factors_per_rule"].append(motif_count)
    if own:
        t["rules_gaining_a_source"] += 1
        t["sources_per_rule"].append(len(own))
    t["candidate_pairs_in_own_cell"] += len(own)
    t["candidate_pairs_with_motif"] += sum(1 for c in own if c["motif"])
    for c in own:
        t["direction"][c["direction"]] += 1
        t["factors"][c["factor"]] += 1
    if any(c["direction"] != "unknown" for c in own):
        t["rules_with_a_directed_candidate"] += 1


def _median(xs: list[int]) -> int | None:
    s = sorted(xs)
    return s[len(s) // 2] if s else None


def finish_tally(t: dict[str, Any], top: int = 10) -> dict[str, Any]:
    return {
        "rules": t["rules"],
        "rules_gaining_a_source": t["rules_gaining_a_source"],
        "by_kind": dict(t["by_kind"]),
        "candidate_pairs_in_own_cell": t["candidate_pairs_in_own_cell"],
        "candidate_pairs_with_motif": t["candidate_pairs_with_motif"],
        "direction_of_candidate_pairs": dict(sorted(t["direction"].items())),
        "rules_with_a_directed_candidate": t["rules_with_a_directed_candidate"],
        "distinct_factors": len(t["factors"]),
        "top_factors_by_rules": t["factors"].most_common(top),
        "sources_per_rule_median": _median(t["sources_per_rule"]),
        "sources_per_rule_max": max(t["sources_per_rule"]) if t["sources_per_rule"] else None,
        "motif_matched_factors_per_rule_median": _median(t["motif_factors_per_rule"]),
        "rules_with_any_motif_match": sum(1 for x in t["motif_factors_per_rule"] if x > 0),
    }


# --- what compile.py must change (not done here: compile.py is another lane's file) ---------------------

COMPILE_SPEC = (
    "compile.py would read the git-ignored assignment table this lane writes "
    "(data/cache/regulators/assignment.tsv.gz: one row per compiled rule and factor occupying its element "
    "in the rule's own cell) and, per program: (1) declare one gene stub per factor it uses, as it "
    "declares target genes, so the factor is an entity a rule can name; (2) for each row whose direction "
    "is activates or inhibits, emit `rule <FACTOR> <direction> <GENE> { through: <element>; when: "
    'cell_type = <cell>; evidence: measured "ENCODE TF ChIP-seq peak cluster, <cell>" occupancy, not '
    "function; direction: curated GO term, factor-level }`, one rule per element, as v0.5 requires; "
    "(3) state no `strength:` it does not have. Two engine facts block a faithful emission today, and "
    "neither is compile.py's to change: the parser reads `strength` only as a number and a rule with no "
    "strength runs at 1.0 (genomeos/ir/model.py), so an unmeasured factor effect would run at full "
    "strength; and a row whose direction is unknown has no grammar token (Albert's decision queue, "
    "'a grammar token for an unresolved direction'), so it must not be emitted as activates or inhibits. "
    "Until both are decided, the honest emission is the evidence as a comment, or a `binds` rule, which "
    "takes no `through` and runs nothing."
)


def biolang_rule(factor: str, rule: CompiledRule, direction: str, basis: str) -> str | None:
    """The v0.5 text one directed row would become, for COMPILE_SPEC's example; None when undirected."""
    if direction not in ("activates", "inhibits"):
        return None
    return (
        f"rule {factor} {direction} {rule.gene} {{ through: {rule.element}; when: cell_type = {rule.cell}; "
        f'evidence: measured "ENCODE TF ChIP-seq peak cluster, {rule.cell}" occupancy, not function; '
        f"direction {basis} }}"
    )


# --- the registered test -------------------------------------------------------------------------------
# Everything below is pure. Constants are module constants, and no function takes a threshold as an
# argument, so the run applies each rule once and has nothing to tune.

TEST_CELLS = ("K562", "HepG2")
KD_ASSAYS = ("CRISPRi RNA-seq", "siRNA RNA-seq", "shRNA RNA-seq", "CRISPR RNA-seq")
TARGET_CLASS = "transcription factor"  # ENCODE target.investigated_as
KD_MAX_RATIO = 0.70  # eligible at >= 30% lower mean TPM of the factor's own gene, knockdown over control
KD_MIN_CONTROL_TPM = 1.0  # the knockdown is assessable only if the factor's gene has >= 1 TPM in control
EXPRESSED_TPM = 1.0  # a universe gene has >= 1 mean TPM in the control replicates
PSEUDO_TPM = 1.0
NEIGHBOUR_BP = 10_000  # genes whose TSS lies this close to the factor's own TSS are left out
PEAK_WINDOW = 50_000  # the baseline counts the factor's peak midpoints within this distance of a TSS
MIN_ARM = 10  # an experiment needs >= 10 genes in the rule arm and >= 10 outside it
STRATA = 10  # deciles of control expression inside the universe
FLOOR_UNITS = 20  # factor-cell units needed for any estimate; a floor, not shown power
BOOTSTRAP_B = 10_000
BOOTSTRAP_SEED = 20261010
PERMUTATIONS = 1_000
PERMUTATION_SEED = 20261011
CONTROL_ALPHA = 0.05

TEST_CONSTANTS = {
    "test_cells": list(TEST_CELLS),
    "kd_assays": list(KD_ASSAYS),
    "target_class": TARGET_CLASS,
    "kd_max_ratio": KD_MAX_RATIO,
    "kd_min_control_tpm": KD_MIN_CONTROL_TPM,
    "expressed_tpm": EXPRESSED_TPM,
    "pseudo_tpm": PSEUDO_TPM,
    "neighbour_bp": NEIGHBOUR_BP,
    "peak_window": PEAK_WINDOW,
    "min_arm": MIN_ARM,
    "strata": STRATA,
    "floor_units": FLOOR_UNITS,
    "bootstrap_b": BOOTSTRAP_B,
    "bootstrap_seed": BOOTSTRAP_SEED,
    "permutations": PERMUTATIONS,
    "permutation_seed": PERMUTATION_SEED,
    "control_alpha": CONTROL_ALPHA,
}


def slim_experiment(exp: dict[str, Any]) -> dict[str, Any]:
    """The fields of an ENCODE experiment object the test reads: identity, target, controls, files."""
    target = exp.get("target") or {}
    controls = [
        c if isinstance(c, str) else c.get("accession", c.get("@id", ""))
        for c in exp.get("possible_controls", [])
    ]
    keys = (
        "accession",
        "output_type",
        "assembly",
        "genome_annotation",
        "biological_replicates",
        "md5sum",
        "file_size",
        "status",
        "file_format",
        "href",
        "date_created",
    )
    files = []
    for f in exp.get("files", []):
        if not isinstance(f, dict):
            continue
        row = {k: f.get(k) for k in keys}
        row["analyses"] = [a.get("@id") if isinstance(a, dict) else a for a in f.get("analyses", [])]
        files.append(row)
    return {
        "accession": exp["accession"],
        "assay_title": exp.get("assay_title"),
        "cell": (exp.get("biosample_ontology") or {}).get("term_name")
        if isinstance(exp.get("biosample_ontology"), dict)
        else None,
        "target": target.get("label") if isinstance(target, dict) else None,
        "investigated_as": target.get("investigated_as", []) if isinstance(target, dict) else [],
        "control_type": exp.get("control_type"),
        "possible_controls": [c.strip("/").split("/")[-1] for c in controls if c],
        "lab": (exp.get("lab") or {}).get("title") if isinstance(exp.get("lab"), dict) else exp.get("lab"),
        "date_released": exp.get("date_released"),
        "status": exp.get("status"),
        "default_analysis": exp.get("default_analysis"),
        "files": files,
    }


def _annotation_rank(annotation: str | None) -> int:
    m = re.fullmatch(r"V(\d+)", annotation or "")
    return int(m.group(1)) if m else -1


def quantifications(exp: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    """Annotation -> one released GRCh38 gene-quantification TSV per biological replicate, preferring the
    file in the experiment's default analysis and then the most recent."""
    groups: dict[tuple[str, tuple], list[dict[str, Any]]] = defaultdict(list)
    for f in exp["files"]:
        if (
            f["output_type"] == "gene quantifications"
            and f["assembly"] == "GRCh38"
            and f["status"] == "released"
            and f["file_format"] == "tsv"
            and _annotation_rank(f["genome_annotation"]) > 0
        ):
            groups[(f["genome_annotation"], tuple(f["biological_replicates"] or ()))].append(f)
    out: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for (ann, _rep), fs in sorted(groups.items()):
        best = max(fs, key=lambda f: (exp.get("default_analysis") in f["analyses"], f["date_created"] or ""))
        out[ann].append(best)
    return dict(out)


def choose_files(kd: dict[str, Any], controls: Sequence[dict[str, Any]]) -> dict[str, Any] | None:
    """The knockdown's and its controls' replicate files on the newest annotation they all share; None
    when there is no control or no shared annotation."""
    if not controls:
        return None
    sets = [quantifications(kd), *(quantifications(c) for c in controls)]
    common = set(sets[0]).intersection(*sets[1:])
    if not common:
        return None
    ann = max(common, key=_annotation_rank)
    return {
        "annotation": ann,
        "knockdown": sets[0][ann],
        "control": [f for s in sets[1:] for f in s[ann]],
    }


def parse_quantification(lines: Iterable[str]) -> dict[str, float]:
    """Unversioned Ensembl gene id -> TPM from an ENCODE RSEM gene-quantification TSV; other rows (spike-
    ins, PAR_Y duplicates) are skipped."""
    out: dict[str, float] = {}
    it = iter(lines)
    header = next(it).rstrip("\n").split("\t")
    gi, ti = header.index("gene_id"), header.index("TPM")
    for line in it:
        f = line.rstrip("\n").split("\t")
        gid = f[gi]
        if not gid.startswith("ENSG") or gid.endswith("_PAR_Y"):
            continue
        out[gid.split(".")[0]] = float(f[ti])
    return out


def mean_tpm(replicates: Sequence[dict[str, float]]) -> dict[str, float]:
    """Gene -> mean TPM over replicates; a gene missing from a replicate counts 0 there."""
    genes: set[str] = set().union(*replicates) if replicates else set()
    n = len(replicates)
    return {g: sum(r.get(g, 0.0) for r in replicates) / n for g in genes}


def knockdown_ratio(factor_gene: str, kd: dict[str, float], ctrl: dict[str, float]) -> float | None:
    """The factor's own mean TPM in knockdown over control; None when not assessable."""
    c = ctrl.get(factor_gene, 0.0)
    if c < KD_MIN_CONTROL_TPM:
        return None
    return kd.get(factor_gene, 0.0) / c


def response(kd: dict[str, float], ctrl: dict[str, float], genes: Iterable[str]) -> dict[str, float]:
    """|log2((TPM_kd + 1) / (TPM_control + 1))| per gene."""
    return {g: abs(math.log2((kd.get(g, 0.0) + PSEUDO_TPM) / (ctrl.get(g, 0.0) + PSEUDO_TPM))) for g in genes}


def deciles(values: Sequence[float]) -> list[int]:
    """Stratum 0..STRATA-1 of each value by rank (tied values share the stratum of their first rank)."""
    n = len(values)
    order = sorted(range(n), key=lambda i: values[i])
    out = [0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and values[order[j + 1]] == values[order[i]]:
            j += 1
        s = min(STRATA - 1, i * STRATA // n)
        for q in range(i, j + 1):
            out[order[q]] = s
        i = j + 1
    return out


def auroc(scores: Sequence[float], labels: Sequence[bool]) -> float | None:
    """P(a labelled score exceeds an unlabelled one), ties counting half; None when a class is empty."""
    pos = sum(1 for x in labels if x)
    neg = len(labels) - pos
    if pos == 0 or neg == 0:
        return None
    order = sorted(range(len(scores)), key=lambda i: scores[i])
    ranks = [0.0] * len(scores)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and scores[order[j + 1]] == scores[order[i]]:
            j += 1
        for q in range(i, j + 1):
            ranks[order[q]] = (i + j) / 2 + 1
        i = j + 1
    rank_sum = sum(r for r, x in zip(ranks, labels, strict=True) if x)
    return (rank_sum - pos * (pos + 1) / 2) / (pos * neg)


def stratified_auroc(scores: Sequence[float], labels: Sequence[bool], strata: Sequence[int]) -> float | None:
    """AUROC over within-stratum pairs only: each stratum's AUROC weighted by its labelled x unlabelled
    pairs. None when no stratum holds both classes."""
    groups: dict[int, list[int]] = defaultdict(list)
    for j, s in enumerate(strata):
        groups[s].append(j)
    num = den = 0.0
    for s in sorted(groups):
        idx = groups[s]
        lab = [labels[j] for j in idx]
        p = sum(lab)
        pairs = p * (len(lab) - p)
        if pairs == 0:
            continue
        a = auroc([scores[j] for j in idx], lab)
        num += a * pairs
        den += pairs
    return num / den if den else None


def near_peaks(
    genes: Sequence[str], tss: dict[str, tuple[str, int]], peaks: dict[str, list[int]]
) -> dict[str, tuple[int, int]]:
    """Gene -> (the factor's peak midpoints within PEAK_WINDOW of its TSS, distance to the nearest
    midpoint on its chromosome, 2**40 when the chromosome has none)."""
    out = {}
    for g in genes:
        chrom, pos = tss[g]
        mids = peaks.get(chrom, [])
        lo = bisect.bisect_left(mids, pos - PEAK_WINDOW)
        hi = bisect.bisect_right(mids, pos + PEAK_WINDOW)
        k = bisect.bisect_left(mids, pos)
        near = [abs(mids[q] - pos) for q in (k - 1, k) if 0 <= q < len(mids)]
        out[g] = (hi - lo, min(near) if near else 1 << 40)
    return out


def baseline_arm(universe: Sequence[str], size: int, near: dict[str, tuple[int, int]]) -> set[str]:
    """The `size` genes of the universe nearest the factor's peaks, ignoring every compiled rule: most
    peaks within PEAK_WINDOW first, then the nearest peak, then the gene id."""
    ranked = sorted(universe, key=lambda g: (-near[g][0], near[g][1], g))
    return set(ranked[:size])


def universe_of(
    rule_genes: Iterable[str],
    ctrl: dict[str, float],
    factor_gene: str,
    tss: dict[str, tuple[str, int]],
) -> list[str]:
    """Genes of the cell's compiled programme the experiment can read: >= EXPRESSED_TPM in control, with a
    TSS, neither the factor's own gene nor within NEIGHBOUR_BP of its TSS. Sorted."""
    fchrom, fpos = tss.get(factor_gene, ("", -(1 << 40)))
    out = []
    for g in sorted(set(rule_genes)):
        if g == factor_gene or g not in tss or ctrl.get(g, 0.0) < EXPRESSED_TPM:
            continue
        chrom, pos = tss[g]
        if chrom == fchrom and abs(pos - fpos) <= NEIGHBOUR_BP:
            continue
        out.append(g)
    return out


def experiment_view(
    universe: Sequence[str], ctrl: dict[str, float], resp: dict[str, float], arm: set[str], base: set[str]
) -> dict[str, Any]:
    """One experiment as the analysis reads it: universe, response, control-expression strata, arms."""
    return {
        "universe": list(universe),
        "y": [resp[g] for g in universe],
        "strata": deciles([math.log2(ctrl[g] + PSEUDO_TPM) for g in universe]),
        "arm": set(arm),
        "base": set(base),
    }


def bootstrap_ci(values: Sequence[float]) -> tuple[float, float]:
    """95% percentile interval of the mean over BOOTSTRAP_B resamples of the units, seeded."""
    rng = random.Random(BOOTSTRAP_SEED)
    n = len(values)
    means = sorted(sum(values[rng.randrange(n)] for _ in range(n)) / n for _ in range(BOOTSTRAP_B))
    return means[int(0.025 * BOOTSTRAP_B)], means[int(0.975 * BOOTSTRAP_B) - 1]


def derangement(items: Sequence[str], rng: random.Random) -> list[str]:
    """A uniformly random permutation of `items` with no fixed point (rejection sampling)."""
    if len(items) < 2:
        raise ValueError("a derangement needs at least two items")
    while True:
        p = list(items)
        rng.shuffle(p)
        if all(a != b for a, b in zip(p, items, strict=True)):
            return p


def decide(n_units: int, d: float | None, lo: float | None, hi: float | None, p: float | None) -> str:
    if n_units < FLOOR_UNITS or d is None or lo is None or hi is None:
        return "stopped_below_floor"
    if p is None or p >= CONTROL_ALPHA:
        return "not_factor_specific"
    if lo > 0:
        return "rules_beat_peaks"
    if hi < 0:
        return "peaks_beat_rules"
    return "no_difference_detected"


OUTCOMES = {
    "stopped_below_floor": (
        "stopped before any comparison: {n} factor-cell units met the registered conditions, below the "
        "floor of {floor}; no estimate is reported and no rule is relaxed"
    ),
    "not_factor_specific": (
        "the compiled regulator assignment did not anticipate a factor's knockdown better than another "
        "factor's assignment did: mean stratified AUROC {a:.4f} for the factor's own assigned genes "
        "against a median of {c:.4f} for another factor's, control p = {p:.4g} over {n} factor-cell units; "
        "the paired difference against genes near the factor's peaks, {d:+.4f} [95% CI {lo:+.4f}, "
        "{hi:+.4f}], is reported and decides nothing"
    ),
    "rules_beat_peaks": (
        "the genes the compiled rules route through a factor's occupied elements changed more on that "
        "factor's knockdown than genes near the factor's peaks: mean paired AUROC difference {d:+.4f} "
        "[95% CI {lo:+.4f}, {hi:+.4f}] over {n} factor-cell units, and more than another factor's "
        "assigned genes (control p = {p:.4g})"
    ),
    "no_difference_detected": (
        "no difference detected: the genes the compiled rules route through a factor's occupied elements "
        "did not change on that factor's knockdown by a detectably different amount from genes near the "
        "factor's peaks (mean paired AUROC difference {d:+.4f} [95% CI {lo:+.4f}, {hi:+.4f}] over {n} "
        "factor-cell units), though they changed more than another factor's assigned genes (control "
        "p = {p:.4g})"
    ),
    "peaks_beat_rules": (
        "genes near a factor's peaks changed more on that factor's knockdown than the genes the compiled "
        "rules route through its occupied elements: mean paired AUROC difference {d:+.4f} [95% CI "
        "{lo:+.4f}, {hi:+.4f}] over {n} factor-cell units (control p = {p:.4g})"
    ),
}


def reading(outcome: str, **numbers: Any) -> str:
    return OUTCOMES[outcome].format(floor=FLOOR_UNITS, **numbers)


def unit_auroc(unit: dict[str, Any], arm_of: str | None = None, base: bool = False) -> float | None:
    """A unit's stratified AUROC, the mean over its experiments: of its own rule arm, of its baseline arm
    (`base`), or of another factor's rule arm in the same cell (`arm_of`), labelled inside this unit's
    universe."""
    vals = []
    for e in unit["experiments"]:
        if base:
            chosen = e["base"]
        elif arm_of is None:
            chosen = e["arm"]
        else:
            chosen = unit["arms"][arm_of]
        a = stratified_auroc(e["y"], [g in chosen for g in e["universe"]], e["strata"])
        if a is not None:
            vals.append(a)
    return sum(vals) / len(vals) if vals else None


def analyse(units: list[dict[str, Any]]) -> dict[str, Any]:
    """The registered analysis over the factor-cell units that met every condition.

    A unit: {factor, cell, experiments: [experiment_view, ...], arms: {factor: genes}}, where `arms`
    holds every unit-factor's rule arm in that cell, for the matched control. The control relabels each
    unit with another unit-factor's rule arm from the same cell (a derangement within each cell type);
    an arm with no gene, or every gene, of the unit's universe scores 0.5 there. A cell type with a
    single unit cannot be deranged within itself, so its unit is set aside before anything is computed
    and reported as `alone_in_cell`."""
    per_cell = Counter(u["cell"] for u in units)
    alone = [{"factor": u["factor"], "cell": u["cell"]} for u in units if per_cell[u["cell"]] < 2]
    units = [u for u in units if per_cell[u["cell"]] >= 2]
    n = len(units)
    if n < FLOOR_UNITS:
        return {
            "outcome": "stopped_below_floor",
            "reading": reading("stopped_below_floor", n=n),
            "units": n,
            "alone_in_cell": alone,
        }
    rules = [unit_auroc(u) for u in units]
    base = [unit_auroc(u, base=True) for u in units]
    if any(v is None for v in (*rules, *base)):
        raise ValueError("a unit with an undefined AUROC was admitted; the admission rule is broken")
    diffs = [a - b for a, b in zip(rules, base, strict=True)]
    d = sum(diffs) / n
    lo, hi = bootstrap_ci(diffs)
    observed = sum(rules) / n
    rng = random.Random(PERMUTATION_SEED)
    by_cell: dict[str, list[int]] = defaultdict(list)
    for i, u in enumerate(units):
        by_cell[u["cell"]].append(i)
    null = []
    for _ in range(PERMUTATIONS):
        vals = [0.5] * n
        for _cell, idx in sorted(by_cell.items()):
            if len(idx) < 2:
                continue
            names = [units[i]["factor"] for i in idx]
            for i, other in zip(idx, derangement(names, rng), strict=True):
                v = unit_auroc(units[i], arm_of=other)
                vals[i] = v if v is not None else 0.5
        null.append(sum(vals) / n)
    p = (1 + sum(1 for x in null if x >= observed)) / (1 + PERMUTATIONS)
    control = sorted(null)[len(null) // 2]
    outcome = decide(n, d, lo, hi, p)
    return {
        "outcome": outcome,
        "reading": reading(outcome, n=n, d=d, lo=lo, hi=hi, p=p, a=observed, c=control),
        "units": n,
        "mean_paired_difference": d,
        "ci95": [lo, hi],
        "mean_auroc_rules": observed,
        "mean_auroc_baseline": sum(base) / n,
        "control_median_mean_auroc": control,
        "control_p": p,
        "alone_in_cell": alone,
        "per_unit": [
            {"factor": u["factor"], "cell": u["cell"], "auroc_rules": a, "auroc_baseline": b}
            for u, a, b in zip(units, rules, base, strict=True)
        ],
    }
