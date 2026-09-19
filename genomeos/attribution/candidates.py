# SPDX-License-Identifier: AGPL-3.0-or-later
"""The syntax candidates read one by one: what each constrained, unannotated block most likely is.

The organiser leaves 69 blocks of the UNKNOWN space constrained on both axes, across 241 mammals
(Zoonomia phyloP) and among 76,156 people (gnomAD Gnocchi), not copies, with no GENCODE gene over
them. This module reads each of them through every layer the repository has already computed or
holds locally, and adds nothing from a model:

- context: the genes on either side (any GENCODE v50 type), which end of each faces the block,
  the nearest coding TSS, the CTCF node and its inferred target (the nearest coding TSS inside
  the node), and the targets the deletion runs named for elements of the same node;
- the human axis's footing: Gnocchi scores whole kilobases, so the kilobases a block touches may
  hold a neighbour's exons, and that constraint is then borrowed, not the block's own;
- constraint inside the block: the 100-vertebrate conserved elements merged into segments, and
  the repeats they sit in (a base in a primate-specific repeat cannot be held across mammals, a
  base in a structured-RNA copy is held because its paralogues are);
- the registry (ENCODE cCREs) and the reader (ENCODE DNase in eleven cell types) on the element
  itself, next to the node's openness in the same cells;
- coding: each conserved segment scored for a stop-free frame with codon usage, held against the
  chromosome's own canonical coding exons, and the segment parser run over the block;
- measured ground truth already local: VISTA, lentiMPRA, ClinVar pathogenic, GWAS lead variants
  and GTEx eQTLs (the last two distilled only near scored elements, so an absence says little).

Every layer is also read on control windows of the same length in nearby intergenic space (no
gene body, no assembly gap), so a candidate's openness or coding score is read against what the
same neighbourhood gives by chance. The reading of each block is a best guess with its basis and
a confidence, never a verdict; a layer that contradicts a class rules the class out.
"""

from __future__ import annotations

import bisect
import gzip
import json
import random
import re
import time
from collections import Counter
from pathlib import Path
from typing import Any

from genomeos.results import RESULTS_DIR, load_result, save_result

CHROMS = tuple([f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"])
CASE = "syntax"
REFERENCE = Path("data/reference")
KNOWLEDGE = Path("data/knowledge")
SIGNALS = RESULTS_DIR / "signals_chr21.json"

MERGE_GAP = 20  # conserved elements this close form one segment
CODING_MIN_SEGMENT = 45  # bp: shorter segments are not scored for coding (15 codons)
STOP_FREE_COVERAGE = 0.9  # a frame counts as open when a stop-free run covers this much of the segment
PARSER_FLANK = 10_000  # sequence either side of the block given to the segment parser
PARSER_CONTROL_EXONS = 40  # internal coding exons per chromosome the parser is asked to find, same flank
CONTROLS_PER_SIDE = 20  # control windows each side of a block
CONTROL_STEP_MIN = 10_000  # control windows start this far apart (or twice the block length)
CONTROL_MAX_TRIES = 200
EXON_CONTROL_SAMPLE = 300  # canonical internal coding exons per chromosome for the coding control
SEED = 20260913
READER_MIN_CELLS = 2  # the reader supports an element open in at least this many cell types ...
READER_MIN_RATIO = 2.0  # ... and at least this many times the control windows' mean
PROMOTER_CLASSES = ("PLS", "DNase-H3K4me3")
ENHANCER_CLASSES = ("dELS", "pELS")
INSULATOR_CLASSES = ("CTCF-only",)
GENE_END_REACH = 3_000  # a conserved segment this close to a gene's 3' end can be its extension
CONFIDENCE_CAP = 0.7  # nothing here is measured on the block itself as a function
REGISTRY_WEIGHT = 0.05  # a cCRE sits in most same-length windows nearby, so the registry adds little alone
READER_WEIGHT = 0.15  # openness above the control windows
MEASURED_WEIGHT = 0.15  # VISTA positive or lentiMPRA active

# repeats whose bases cannot be held across placental mammals: they are younger than the clade
PRIMATE_SPECIFIC = re.compile(r"^(Alu|FLAM|FRAM|L1PA|L1PB|L1P\d|SVA|HERV|LTR12|LTR5|MER4\d|THE1|PABL)")
STRUCTURED_RNA_FAMILIES = {"srpRNA", "tRNA", "snRNA", "rRNA", "scRNA", "7SK", "RNA"}

EVIDENCE = {
    "context": "curated: GENCODE v50 genes (every type); inferred: CTCF-only nodes, nearest coding TSS",
    "node_targets": "predicted: AlphaGenome deletion runs already scored (no new calls)",
    "registry": "curated: ENCODE SCREEN cCREs v3",
    "reader": "experimental: ENCODE DNase-seq narrowPeak, eleven cell types; control windows alongside",
    "constraint": "experimental-derived: UCSC phastConsElements100way (local cache)",
    "repeats": "curated: UCSC RepeatMasker",
    "coding": (
        "predicted: codon log-odds of the best stop-free frame, against the chromosome's canonical "
        "coding exons and control segments; segment parser v1 (learned PWMs, Viterbi)"
    ),
    "ground_truth": (
        "experimental: VISTA e11.5, ENCODE4 lentiMPRA, GTEx v8 eQTLs; curated: ClinVar pathogenic, "
        "GWAS Catalog lead variants"
    ),
    "reading": "inferred: rules over the layers above, a best guess with a confidence",
}

CLASSES = (
    "alignment_artefact",
    "structured_rna_copy",
    "coding_copy",
    "promoter",
    "regulatory",
    "coding_exon",
    "transcript_extension",
    "unexplained",
)
GROUP = {
    "alignment_artefact": "artefact",
    "structured_rna_copy": "rna",
    "coding_copy": "coding",
    "promoter": "regulatory",
    "regulatory": "regulatory",
    "coding_exon": "coding",
    "transcript_extension": "rna",
    "unexplained": "unexplained",
}


# ---- intervals ---------------------------------------------------------------------------


class Intervals:
    """Sorted half-open intervals with a payload; overlap queries by bisection."""

    def __init__(self, items: list[tuple[int, int, Any]]) -> None:
        self.items = sorted(items, key=lambda x: (x[0], x[1]))
        self.starts = [x[0] for x in self.items]
        self.max_len = max((e - s for s, e, _ in self.items), default=0)

    def __len__(self) -> int:
        return len(self.items)

    def overlapping(self, start: int, end: int) -> list[tuple[int, int, Any]]:
        i = bisect.bisect_left(self.starts, start - self.max_len)
        out = []
        for s, e, v in self.items[i:]:
            if s >= end:
                break
            if e > start:
                out.append((s, e, v))
        return out

    def any(self, start: int, end: int) -> bool:
        return bool(self.overlapping(start, end))

    def covered(self, start: int, end: int) -> int:
        return sum(
            b - a
            for a, b in merge([(max(s, start), min(e, end)) for s, e, _ in self.overlapping(start, end)])
        )


def merge(intervals: list[tuple[int, int]], gap: int = 0) -> list[tuple[int, int]]:
    out: list[list[int]] = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1] + gap:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def overlap_bp(a: tuple[int, int], b: tuple[int, int]) -> int:
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


# ---- coding ------------------------------------------------------------------------------

STOPS = {"TAA", "TAG", "TGA"}
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


_BASES = "TCAG"
_AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {
    a + b + c: _AMINO[16 * i + 4 * j + k]
    for i, a in enumerate(_BASES)
    for j, b in enumerate(_BASES)
    for k, c in enumerate(_BASES)
}
KMER = 8  # peptide k-mer shared with a known protein: 20^8 makes a chance match rare
PARALOGUE_MIN_SHARE = 0.3  # share of a query's k-mers found in one protein that reads as a copy of it
PARALOGUE_MIN_KMERS = 20  # a shorter peptide cannot carry a copy call
KMER_MIN_DISTINCT = (
    5  # distinct residues a k-mer needs: collagen-like and other low-complexity runs match anything
)
CODING_TAIL_P = 0.05  # coding-like segments beyond what the control rate gives by chance (binomial tail)


def translate(dna: str) -> str:
    return "".join(CODE.get(dna[i : i + 3], "X") for i in range(0, len(dna) - 2, 3))


def frame_score(seq: str, codon_lo: dict[str, float]) -> dict[str, Any]:
    """The segment's best open frame over six: the longest stop-free run of codons in each frame, and
    its mean codon log-odds where that run covers STOP_FREE_COVERAGE of the segment. A segment with
    no open frame cannot be a coding exon read through, whatever its conservation."""
    s = seq.upper()
    n = len(s)
    best_cov, best_score, best_peptide = 0.0, None, ""
    for strand_seq in (s, revcomp(s)):
        for f in range(3):
            run, run_start, longest = 0, f, (0, f)
            for i in range(f, n - 2, 3):
                if strand_seq[i : i + 3] in STOPS:
                    run = 0
                    run_start = i + 3
                    continue
                run += 1
                if run > longest[0]:
                    longest = (run, run_start)
            cov = longest[0] * 3 / n if n else 0.0
            best_cov = max(best_cov, cov)
            if cov >= STOP_FREE_COVERAGE and longest[0]:
                a = longest[1]
                codons = [strand_seq[a + 3 * k : a + 3 * k + 3] for k in range(longest[0])]
                sc = sum(codon_lo.get(c, 0.0) for c in codons) / len(codons)
                if best_score is None or sc > best_score:
                    best_score = sc
                    best_peptide = "".join(CODE.get(c, "X") for c in codons)
    return {
        "open_frame": best_score is not None,
        "score": round(best_score, 4) if best_score is not None else None,
        "coverage": round(best_cov, 3),
        "peptide": best_peptide,
    }


def paralogues(
    queries: dict[str, str], proteins: Path = KNOWLEDGE / "proteins", k: int = KMER
) -> dict[str, dict[str, Any]]:
    """For each query peptide, the known human protein (UniProt, local) sharing most of its k-mers.
    One pass over the proteome, checking membership in the queries' k-mer set only."""
    kmers: dict[str, set[str]] = {}
    for key, pep in queries.items():
        kmers[key] = {
            pep[i : i + k]
            for i in range(len(pep) - k + 1)
            if not set(pep[i : i + k]) & {"X", "*"} and len(set(pep[i : i + k])) >= KMER_MIN_DISTINCT
        }
    wanted: set[str] = set().union(*kmers.values()) if kmers else set()
    if not wanted or not proteins.exists():
        return {}
    hits: dict[str, list[tuple[str, str]]] = {}
    for p in sorted(proteins.glob("*.json")):
        try:
            d = json.loads(p.read_text())
        except (OSError, ValueError):
            continue
        items = ((d.get("sections") or {}).get("identity") or {}).get("items") or {}
        seq = items.get("sequence") or ""
        if not seq:
            continue
        found = {seq[i : i + k] for i in range(len(seq) - k + 1)} & wanted
        for km in found:
            hits.setdefault(km, []).append((d.get("gene") or p.stem, items.get("accession") or ""))
    out: dict[str, dict[str, Any]] = {}
    for key, ks in kmers.items():
        per: Counter = Counter()
        for km in ks:
            for gene_acc in set(hits.get(km, [])):
                per[gene_acc] += 1
        if per:
            (gene, acc), n = per.most_common(1)[0]
            out[key] = {
                "protein": gene,
                "accession": acc,
                "shared_kmers": n,
                "query_kmers": len(ks),
                "share": round(n / max(1, len(ks)), 3),
            }
        else:
            out[key] = {"protein": None, "shared_kmers": 0, "query_kmers": len(ks), "share": 0.0}
    return out


def percentile(value: float, pool: list[float]) -> float | None:
    if not pool:
        return None
    return round(sum(1 for x in pool if x <= value) / len(pool), 3)


def quantile(pool: list[float], q: float) -> float | None:
    if not pool:
        return None
    xs = sorted(pool)
    return xs[min(len(xs) - 1, max(0, int(round(q * (len(xs) - 1)))))]


# ---- per-chromosome evidence ---------------------------------------------------------------


def candidates(results_dir: Path = RESULTS_DIR, case: str = CASE) -> list[dict[str, Any]]:
    """The organiser's constrained_unknown blocks of one case that are not copies, every chromosome."""
    out = []
    for c in CHROMS:
        r = load_result(f"organised_{c}", results_dir)
        if not r:
            continue
        rows = [x for x in r.get("candidates", []) if x.get("case") == case and not x.get("copy")]
        expected = ((r.get("real_unknown") or {}).get("by_case") or {}).get(case, {}).get("blocks", 0)
        if len(rows) < expected:
            raise ValueError(
                f"organised_{c} lists {len(rows)} {case} candidates of {expected}; "
                "rerun organise with a larger top"
            )
        out.extend({**x, "chrom": c} for x in rows)
    return out


def _bed(path: Path):
    if not path.exists():
        return
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            yield line.rstrip("\n").split("\t")


def _tss(g) -> int:
    return g.locus.end - 1 if g.locus.strand.value == "-" else g.locus.start


def cells_available(results_dir: Path = RESULTS_DIR) -> list[str]:
    return sorted(
        p.name[len("dnase_") : -len("_chr1.bed.gz")] for p in results_dir.glob("dnase_*_chr1.bed.gz")
    )


def deletion_elements(chrom: str, results_dir: Path = RESULTS_DIR) -> list[dict[str, Any]]:
    """Elements the deletion runs have already scored on this chromosome. The whole-chromosome run is
    read only when its summary says it is complete; the samples are read inline."""
    from genomeos.attribution.targets import run_elements

    out: dict[str, dict[str, Any]] = {}
    for name in ("enhancer_targets_all", "constrained_targets", "enhancer_targets"):
        r = load_result(f"{name}_{chrom}", results_dir) or {}
        if name == "enhancer_targets_all" and not r.get("complete"):
            continue
        try:
            rows = run_elements(name, chrom, results_dir)
        except FileNotFoundError:
            continue
        for e in rows:
            out.setdefault(e["id"], {**e, "origin": name})
    return list(out.values())


class Chromosome:
    """Everything local about one chromosome that the candidate reading needs, loaded once."""

    def __init__(self, chrom: str, results_dir: Path = RESULTS_DIR, reference: Path = REFERENCE) -> None:
        from genomeos.genome import Annotation
        from genomeos.genome.index import IndexedGenome

        self.chrom = chrom
        self.results_dir = results_dir
        self.annotation = Annotation.from_gff3(reference / f"gencode_v50_{chrom}.gff3.gz", {chrom})
        self.genes = sorted(self.annotation.genes.values(), key=lambda g: g.locus.start)
        self.gene_ends = sorted((g.locus.end, i) for i, g in enumerate(self.genes))
        self.bodies = merge([(g.locus.start, g.locus.end) for g in self.genes])
        self.body_starts = [s for s, _ in self.bodies]
        exons = []
        for g in self.genes:
            for t in g.transcripts.values():
                exons.extend((e.start, e.end) for e in t.exons)
        self.exons = merge(exons)
        self.exon_starts = [s for s, _ in self.exons]
        self.coding = [g for g in self.genes if g.type == "protein_coding"]
        self.coding_tss = sorted((_tss(g), g.symbol) for g in self.coding)
        self.by_id = {g.id.split(".")[0]: g.symbol for g in self.genes}
        self.genome = IndexedGenome(reference / f"{chrom}.fa")
        self.length = self.genome.lengths[chrom]
        self.ccres = Intervals(
            [(int(f[1]), int(f[2]), (f[3], f[4])) for f in _bed(results_dir / f"ccres_{chrom}.bed.gz")]
        )
        self.conserved = Intervals(
            [
                (int(f[1]), int(f[2]), int(f[3]))
                for f in _bed(KNOWLEDGE / "constraint" / f"phastConsElements100way_{chrom}.bed.gz")
            ]
        )
        self.repeats = Intervals(
            [(int(f[0]), int(f[1]), (f[2], f[3], f[4])) for f in _bed(results_dir / f"rmsk_{chrom}.bed.gz")]
        )
        self.cells = cells_available(results_dir)
        self.dnase = {
            c: Intervals(
                [(int(f[0]), int(f[1]), float(f[2])) for f in _bed(results_dir / f"dnase_{c}_{chrom}.bed.gz")]
            )
            for c in self.cells
        }
        self.domains = sorted(
            (load_result(f"domains_{chrom}", results_dir) or {}).get("domains", []), key=lambda d: d["start"]
        )
        self.domain_starts = [d["start"] for d in self.domains]
        self.node_open: dict[str, dict[str, bool]] = {}
        for c in self.cells:
            r = load_result(f"reader_{c}_{chrom}", results_dir) or {}
            table = r.get("node_table") or []
            dens = sorted(n["peaks_per_100kb"] for n in table)
            median = dens[len(dens) // 2] if dens else 0.0
            self.node_open[c] = {n["id"]: n["peaks_per_100kb"] >= max(1.0, median) for n in table}
        self.deleted = deletion_elements(chrom, results_dir)
        self._codon_lo: dict[str, float] | None = None
        self._exon_scores: list[float] | None = None

    def close(self) -> None:
        self.genome.close()

    # sequence ------------------------------------------------------------------------------
    def seq(self, start: int, end: int) -> str:
        from genomeos.coords import Locus

        start, end = max(0, start), min(self.length, end)
        return str(self.genome.fetch(Locus(self.chrom, start, end))).upper()

    def canonical_cds(self) -> list[list[tuple[int, int]]]:
        out = []
        for g in self.coding:
            for t in g.transcripts.values():
                if t.cds and "Ensembl_canonical" in t.tags:
                    out.append((sorted((c.start, c.end) for c, _ in t.cds), g.locus.strand.value))
        return out

    @property
    def codon_lo(self) -> dict[str, float]:
        """Codon log-odds learned from this chromosome's canonical coding sequence against a 3 Mb sample
        of the chromosome (the segment parser's own table)."""
        if self._codon_lo is None:
            from genomeos.genome.segments import codon_log_odds

            cds = []
            for segs, strand in self.canonical_cds():
                s = "".join(self.seq(a, b) for a, b in segs)
                cds.append(revcomp(s) if strand == "-" else s)
            mid = self.length // 3
            self._codon_lo = codon_log_odds(cds, self.seq(mid, mid + 3_000_000))
        return self._codon_lo

    @property
    def exon_scores(self) -> list[float]:
        """The positive control: canonical internal coding exons, scored blind to strand and frame."""
        if self._exon_scores is None:
            internal = []
            for segs, _ in self.canonical_cds():
                internal.extend((a, b) for a, b in segs[1:-1] if CODING_MIN_SEGMENT <= b - a <= 1_000)
            rng = random.Random(SEED)
            sample = rng.sample(internal, min(EXON_CONTROL_SAMPLE, len(internal)))
            scores = [frame_score(self.seq(a, b), self.codon_lo)["score"] for a, b in sample]
            self._exon_scores = [s for s in scores if s is not None]
        return self._exon_scores

    def parser_sensitivity(self, parser) -> dict[str, Any]:
        """The parser's positive control on this chromosome: internal coding exons it recovers when given
        the same flank as a candidate. Without it, "no exon predicted" would read as evidence."""
        internal = [(a, b) for segs, _ in self.canonical_cds() for a, b in segs[1:-1]]
        rng = random.Random(SEED + 1)
        sample = rng.sample(internal, min(PARSER_CONTROL_EXONS, len(internal)))
        found = 0
        for a, b in sample:
            wa, wb = max(0, a - PARSER_FLANK), min(self.length, b + PARSER_FLANK)
            preds = parser.parse(self.chrom, self.seq(wa, wb), wa)
            if any(overlap_bp((wa + x, wa + y), (a, b)) for p in preds for x, y in p.exons):
                found += 1
        return {
            "exons": len(sample),
            "recovered": found,
            "sensitivity": round(found / max(1, len(sample)), 3),
        }

    # lookups -------------------------------------------------------------------------------
    def in_gene_body(self, start: int, end: int) -> bool:
        i = bisect.bisect_right(self.body_starts, end) - 1
        while i >= 0 and self.bodies[i][1] > start:
            if self.bodies[i][0] < end:
                return True
            i -= 1
        return False

    def exon_bp(self, start: int, end: int) -> int:
        i = max(0, bisect.bisect_right(self.exon_starts, start) - 1)
        tot = 0
        for s, e in self.exons[i:]:
            if s >= end:
                break
            tot += overlap_bp((s, e), (start, end))
        return tot

    def node_at(self, pos: int) -> dict[str, Any] | None:
        i = bisect.bisect_right(self.domain_starts, pos) - 1
        if 0 <= i < len(self.domains) and self.domains[i]["start"] <= pos < self.domains[i]["end"]:
            return self.domains[i]
        return None


# ---- one block -----------------------------------------------------------------------------


def flank_genes(ch: Chromosome, start: int, end: int) -> dict[str, Any]:
    """The nearest gene on each side and which of its ends faces the block."""
    left = None
    i = bisect.bisect_right([e for e, _ in ch.gene_ends], start) - 1
    if i >= 0:
        left = ch.genes[ch.gene_ends[i][1]]
    j = bisect.bisect_left([g.locus.start for g in ch.genes], end)
    right = ch.genes[j] if j < len(ch.genes) else None

    def side(g, facing_left: bool) -> dict[str, Any] | None:
        if g is None:
            return None
        plus = g.locus.strand.value != "-"
        # the left gene faces the block with its end; that end is the 3' end on the plus strand
        three_prime = plus if facing_left else not plus
        dist = start - g.locus.end if facing_left else g.locus.start - end
        return {
            "gene": g.symbol,
            "type": g.type,
            "strand": g.locus.strand.value,
            "distance": max(0, dist),
            "facing": "3' end" if three_prime else "5' end (TSS)",
        }

    return {"left": side(left, True), "right": side(right, False)}


def nearest_coding_tss(ch: Chromosome, pos: int, lo: int | None = None, hi: int | None = None):
    best = None
    i = bisect.bisect_left(ch.coding_tss, (pos, ""))
    for j in range(max(0, i - 3), min(len(ch.coding_tss), i + 3)):
        t, sym = ch.coding_tss[j]
        if lo is not None and not (lo <= t < hi):
            continue
        d = abs(t - pos)
        if best is None or d < best[1]:
            best = (sym, d)
    if best is None and lo is not None:
        inside = [(abs(t - pos), sym) for t, sym in ch.coding_tss if lo <= t < hi]
        if inside:
            d, sym = min(inside)
            best = (sym, d)
    return {"gene": best[0], "distance": best[1]} if best else None


def segments(ch: Chromosome, start: int, end: int) -> list[dict[str, Any]]:
    hits = ch.conserved.overlapping(start, end)
    merged = merge([(max(s, start), min(e, end)) for s, e, _ in hits], MERGE_GAP)
    out = []
    for a, b in merged:
        lod = max((v for s, e, v in hits if s < b and e > a), default=0)
        out.append({"start": a, "end": b, "length": b - a, "max_lod": lod})
    return out


def repeat_context(ch: Chromosome, start: int, end: int, segs: list[dict[str, Any]]) -> dict[str, Any]:
    reps = ch.repeats.overlapping(start, end)
    covered = merge([(max(s, start), min(e, end)) for s, e, _ in reps])
    groups: Counter = Counter()
    families: Counter = Counter()
    for seg in segs:
        for p in range(seg["start"], seg["end"]):
            hit = next(((cls, fam, name) for s, e, (cls, fam, name) in reps if s <= p < e), None)
            if hit is None:
                groups["unique"] += 1
                continue
            cls, fam, name = hit
            families[f"{fam}:{name}"] += 1
            if fam in STRUCTURED_RNA_FAMILIES or cls in STRUCTURED_RNA_FAMILIES:
                groups["structured_rna"] += 1
            elif PRIMATE_SPECIFIC.match(name):
                groups["primate_specific"] += 1
            elif cls in ("Simple_repeat", "Low_complexity"):
                groups["simple"] += 1
            else:
                groups["ancient_repeat"] += 1
    cons_bp = sum(s["length"] for s in segs)
    return {
        "repeat_fraction": round(sum(b - a for a, b in covered) / max(1, end - start), 3),
        "conserved_bp": cons_bp,
        "conserved_bp_by_repeat_age": dict(groups),
        "conserved_families": [f for f, _ in families.most_common(4)],
    }


def reader(
    ch: Chromosome, start: int, end: int, segs: list[dict[str, Any]], node_id: str | None
) -> dict[str, Any]:
    open_element, open_conserved, open_node = [], [], []
    for c in ch.cells:
        idx = ch.dnase[c]
        if idx.any(start, end):
            open_element.append(c)
        if any(idx.any(s["start"], s["end"]) for s in segs):
            open_conserved.append(c)
        if node_id and ch.node_open[c].get(node_id):
            open_node.append(c)
    return {
        "cells": len(ch.cells),
        "open_on_element": open_element,
        "open_on_conserved": open_conserved,
        "node_open": open_node,
        "node_open_element_closed": sorted(set(open_node) - set(open_element)),
    }


def parser_hits(ch: Chromosome, start: int, end: int, segs: list[dict[str, Any]], parser) -> dict[str, Any]:
    a, b = max(0, start - PARSER_FLANK), min(ch.length, end + PARSER_FLANK)
    seq = ch.seq(a, b)
    preds = parser.parse(ch.chrom, seq, a)
    exons_in, on_conserved, best = 0, 0, None
    for p in preds:
        ex = [(a + x, a + y) for x, y in p.exons]
        inside = [e for e in ex if overlap_bp(e, (start, end))]
        if not inside:
            continue
        exons_in += len(inside)
        cons = [e for e in inside if any(overlap_bp(e, (s["start"], s["end"])) for s in segs)]
        on_conserved += len(cons)
        if cons and (best is None or p.score > best["score"]):
            outside = [e for e in ex if not overlap_bp(e, (start, end)) and ch.exon_bp(*e) > 0]
            joined = sorted(
                {
                    g.symbol
                    for e in outside
                    for g in ch.genes
                    if g.locus.start < e[1] and g.locus.end > e[0] and g.locus.strand.value == p.strand.value
                }
            )
            best = {
                "strand": p.strand.value,
                "start": a + p.start,
                "end": a + p.end,
                "exons": len(ex),
                "coding_bp": sum(y - x for x, y in ex),
                "score": round(p.score, 1),
                "joins_annotated_exons_of": joined,
                "_exons": [(max(x, start), min(y, end)) for x, y in inside],
            }
    if best is not None:
        # translate only the part inside the block: outside it the structure may be a known gene
        dna = "".join(ch.seq(x, y) for x, y in sorted(best.pop("_exons")))
        best["peptide"] = translate(revcomp(dna) if best["strand"] == "-" else dna)
    return {"exons_in_block": exons_in, "exons_on_conserved": on_conserved, "best_on_conserved": best}


def coding(ch: Chromosome, segs: list[dict[str, Any]]) -> dict[str, Any]:
    scored = []
    for s in segs:
        if s["length"] < CODING_MIN_SEGMENT:
            continue
        fs = frame_score(ch.seq(s["start"], s["end"]), ch.codon_lo)
        reps = [
            (a, b, v)
            for a, b, v in ch.repeats.overlapping(s["start"], s["end"])
            if v[0] not in ("Simple_repeat", "Low_complexity")
        ]
        cov = sum(x - y for y, x in merge([(max(a, s["start"]), min(b, s["end"])) for a, b, _ in reps]))
        top = max(reps, key=lambda r: overlap_bp((r[0], r[1]), (s["start"], s["end"])), default=None)
        in_repeat = f"{top[2][1]}:{top[2][2]}" if top and cov >= 0.5 * s["length"] else None
        scored.append({"start": s["start"], "end": s["end"], **fs, "in_repeat": in_repeat})
    return {"segments_scored": len(scored), "segments": scored}


def ground_truth(ch: Chromosome, start: int, end: int, truth: dict[str, Any]) -> dict[str, Any]:
    vista = [v for v in truth["vista"] if v.start < end and v.end > start]
    mpra = [m for m in truth["mpra"] if m["start"] < end and m["end"] > start]
    clin = [r for r in truth["clinvar"] if start <= r["pos"] - 1 < end]
    gwas = [r for r in truth["gwas"] if start - 1_000 <= r["pos"] - 1 < end + 1_000]
    eqtl = [r for r in truth["eqtl"] if start <= r["pos"] - 1 < end]
    return {
        "vista": [{"id": v.id, "status": v.status, "tissues": list(v.tissues)} for v in vista],
        "mpra": [{"key": m["key"], "activity": m["activity"], "active": m["active"]} for m in mpra],
        "clinvar_pathogenic": [
            {
                "pos": r["pos"],
                "gene": r["gene"],
                "significance": r["significance"],
                "consequence": r["consequence"],
            }
            for r in clin
        ],
        "gwas_within_1kb": sorted({(r["rs"], r["trait"]) for r in gwas}),
        "eqtl_genes": sorted({ch.by_id.get(r["gene_id"].split(".")[0], r["gene_id"]) for r in eqtl}),
        "eqtl_tissues": len({r["tissue"] for r in eqtl}),
    }


def load_truth(chrom: str) -> dict[str, Any]:
    from genomeos.attribution import eqtl, gwas, mpra, vista

    vpath = vista.locus_path()
    vrows = []
    if vpath.exists():
        with gzip.open(vpath, "rt") as fh:
            vrows = vista.parse_loci(fh, chrom)
    clin = []
    for f in _bed(KNOWLEDGE / "clinvar" / "pathogenic.tsv.gz"):
        if f[0] == chrom:
            clin.append(
                {
                    "pos": int(f[1]),
                    "gene": f[4],
                    "significance": f[5].replace("_", " "),
                    "consequence": f[11] if len(f) > 11 else "",
                }
            )
    seen: set = set()
    g_rows = []
    for rows in gwas.load_hits().values():
        for r in rows:
            k = (r["chrom"], r["pos"], r["rs"], r["trait"])
            if r["chrom"] == chrom and k not in seen:
                seen.add(k)
                g_rows.append(r)
    e_rows = []
    seen = set()
    for rows in eqtl.load_hits(chrom=chrom).values():
        for r in rows:
            k = (r["tissue"], r["pos"], r["gene_id"])
            if k not in seen:
                seen.add(k)
                e_rows.append(r)
    return {"vista": vrows, "mpra": mpra.load_rows(chrom), "clinvar": clin, "gwas": g_rows, "eqtl": e_rows}


def node_context(ch: Chromosome, start: int, end: int) -> dict[str, Any]:
    mid = (start + end) // 2
    node = ch.node_at(mid)
    if node is None:
        return {"node": None}
    inferred = nearest_coding_tss(ch, mid, node["start"], node["end"])
    inside = [e for e in ch.deleted if e.get("domain") == node["id"]]
    named = Counter((e.get("predicted_coding") or {}).get("gene") for e in inside)
    named.pop(None, None)
    on_block = [
        {
            "id": e["id"],
            "target": (e.get("predicted_coding") or {}).get("gene"),
            "log2": (e.get("predicted_coding") or {}).get("log2_fold_change"),
        }
        for e in inside
        if e["start"] < end and e["end"] > start
    ]
    return {
        "node": node["id"],
        "node_start": node["start"],
        "node_end": node["end"],
        "node_coding_genes": node.get("coding_genes"),
        "inferred_target": inferred,
        "deleted_elements_in_node": len(inside),
        "node_targets_named": [{"gene": g, "elements": n} for g, n in named.most_common(3)],
        "deleted_elements_on_block": on_block,
    }


def human_axis_footing(ch: Chromosome, start: int, end: int) -> dict[str, Any]:
    """Gnocchi counts whole kilobases: the ones this block touches, and how much of each is a gene's exon."""
    k0, k1 = (start // 1000) * 1000, -(-end // 1000) * 1000
    windows = [(k, k + 1000) for k in range(k0, k1, 1000)]
    with_exon = sum(1 for w in windows if ch.exon_bp(*w) > 0)
    return {
        "kilobases": len(windows),
        "kilobases_with_exon_bases": with_exon,
        "exon_bp_in_kilobases": sum(ch.exon_bp(*w) for w in windows),
    }


def control_windows(ch: Chromosome, start: int, end: int) -> list[tuple[int, int]]:
    """Same-length windows in nearby intergenic space: no gene body, no assembly gap, not the block."""
    length = end - start
    step = max(CONTROL_STEP_MIN, 2 * length)
    out = []
    for sign in (-1, 1):
        kept, k = 0, 0
        while kept < CONTROLS_PER_SIDE and k < CONTROL_MAX_TRIES:
            k += 1
            a = start + sign * k * step
            b = a + length
            if a < 0 or b > ch.length or ch.in_gene_body(a, b):
                continue
            s = ch.seq(a, b)
            if s.count("N") > 0.01 * length:
                continue
            out.append((a, b))
            kept += 1
    return out


def read_window(
    ch: Chromosome, start: int, end: int, parser=None, node_id: str | None = None
) -> dict[str, Any]:
    """The layers that controls share with candidates."""
    segs = segments(ch, start, end)
    ccre = [v for _, _, v in ch.ccres.overlapping(start, end)]
    on_cons = {v for s in segs for _, _, v in ch.ccres.overlapping(s["start"], s["end"])}
    cons_bp = sum(s["length"] for s in segs)
    out = {
        "segments": segs,
        "ccre_classes": sorted({cls for _, cls in ccre}),
        "ccres": [i for i, _ in ccre],
        "ccre_classes_on_conserved": sorted({cls for _, cls in on_cons}),
        "reader": reader(ch, start, end, segs, node_id),
        "coding": coding(ch, segs),
        "bp": {
            "window": end - start,
            "in_ccre": ch.ccres.covered(start, end),
            "conserved": cons_bp,
            "conserved_in_ccre": sum(ch.ccres.covered(s["start"], s["end"]) for s in segs),
            # summed over cell types: divide by the cell count for a per-cell fraction
            "open": sum(ch.dnase[c].covered(start, end) for c in ch.cells),
            "conserved_open": sum(ch.dnase[c].covered(s["start"], s["end"]) for c in ch.cells for s in segs),
        },
    }
    if parser is not None:
        out["parser"] = parser_hits(ch, start, end, segs, parser)
    return out


def read_block(ch: Chromosome, row: dict[str, Any], truth: dict[str, Any], parser) -> dict[str, Any]:
    start, end = row["start"], row["end"]
    ctx = node_context(ch, start, end)
    w = read_window(ch, start, end, parser, ctx.get("node"))
    controls = [read_window(ch, a, b, parser) for a, b in control_windows(ch, start, end)]
    exon_pool = ch.exon_scores
    exon_median = quantile(exon_pool, 0.5)
    for s in w["coding"]["segments"]:
        s["exon_percentile"] = percentile(s["score"], exon_pool) if s["score"] is not None else None
        s["coding_like"] = bool(
            s["score"] is not None and exon_median is not None and s["score"] >= exon_median
        )
    n_ctrl = max(1, len(controls))
    ctrl_open = sum(len(c["reader"]["open_on_element"]) for c in controls) / n_ctrl
    ctrl_segments = [s for c in controls for s in c["coding"]["segments"]]
    mid = (start + end) // 2
    return {
        "chrom": ch.chrom,
        "start": start,
        "end": end,
        "length": end - start,
        "class": row.get("class"),
        "mammal_fraction": row.get("mammal_fraction"),
        "human_fraction": row.get("human_fraction"),
        "conserved_elements": row.get("conserved_elements"),
        "flanks": flank_genes(ch, start, end),
        "nearest_coding_tss": nearest_coding_tss(ch, mid),
        **ctx,
        "human_axis_footing": human_axis_footing(ch, start, end),
        "segments": len(w["segments"]),
        "segment_spans": [[s["start"], s["end"]] for s in w["segments"]],
        "longest_segment": max((s["length"] for s in w["segments"]), default=0),
        "repeats": repeat_context(ch, start, end, w["segments"]),
        "ccre_classes": w["ccre_classes"],
        "ccre_classes_on_conserved": w["ccre_classes_on_conserved"],
        "ccres": w["ccres"],
        "reader": w["reader"],
        "coding": w["coding"],
        "parser": w["parser"],
        "bp": w["bp"],
        "ground_truth": ground_truth(ch, start, end, truth),
        "controls": {
            "windows": len(controls),
            "mean_cells_open": round(ctrl_open, 2),
            "with_ccre": round(sum(1 for c in controls if c["ccres"]) / n_ctrl, 3),
            "with_promoter_ccre": round(
                sum(1 for c in controls if set(c["ccre_classes"]) & set(PROMOTER_CLASSES)) / n_ctrl, 3
            ),
            "bp": {k: sum(c["bp"][k] for c in controls) for k in w["bp"]},
            "with_parser_exon_on_conserved": round(
                sum(1 for c in controls if c["parser"]["exons_on_conserved"]) / n_ctrl, 3
            ),
            "segments_scored": len(ctrl_segments),
            "segments_open_frame": sum(1 for x in ctrl_segments if x["score"] is not None),
            "segments_coding_like": sum(
                1
                for x in ctrl_segments
                if x["score"] is not None and exon_median is not None and x["score"] >= exon_median
            ),
            "_per_window": [
                {
                    "cells_open": len(c["reader"]["open_on_element"]),
                    "ccre": bool(c["ccres"]),
                    "parser_exon_on_conserved": bool(c["parser"]["exons_on_conserved"]),
                    "segments": len(c["segments"]),
                }
                for c in controls
            ],
        },
        "exon_control": {"exons": len(exon_pool), "median_score": exon_median},
    }


# ---- the reading ---------------------------------------------------------------------------


def _gene_end_segment(ev: dict[str, Any]) -> dict[str, Any] | None:
    """A flanking gene whose 3' end faces the block with a conserved segment within reach of it."""
    for side in ("left", "right"):
        g = ev["flanks"].get(side)
        if not g or g["facing"] != "3' end" or g["type"] != "protein_coding":
            continue
        edge = ev["start"] - g["distance"] if side == "left" else ev["end"] + g["distance"]
        if any(min(abs(a - edge), abs(b - edge)) <= GENE_END_REACH for a, b in ev["segment_spans"]):
            return g
    return None


def reading(ev: dict[str, Any]) -> dict[str, Any]:
    """A best guess for one block from its layers: class, confidence, basis, and what was ruled out."""
    basis: list[str] = []
    notes: list[str] = []
    ruled_out: list[str] = []
    cells = ev["reader"]["open_on_element"]
    expected = ev["controls"]["mean_cells_open"]
    reader_ok = len(cells) >= READER_MIN_CELLS and len(cells) >= READER_MIN_RATIO * max(expected, 0.5)
    # the registry counts where it sits on the constrained bases, not merely somewhere in the block:
    # a cCRE lies in most same-length windows of these neighbourhoods
    classes = set(ev["ccre_classes_on_conserved"])
    if set(ev["ccre_classes"]) - classes:
        notes.append(
            f"registry: {'/'.join(sorted(set(ev['ccre_classes']) - classes))} "
            "in the block but off its conserved bases"
        )
    promoter = bool(classes & set(PROMOTER_CLASSES))
    enhancer = bool(classes & set(ENHANCER_CLASSES))
    insulator = bool(classes & set(INSULATOR_CLASSES))
    gt = ev["ground_truth"]
    vista_pos = [v["id"] for v in gt["vista"] if v["status"] == "positive"]
    mpra_active = [m["key"] for m in gt["mpra"] if any(m["active"].values())]
    ages = ev["repeats"]["conserved_bp_by_repeat_age"]
    cons_bp = max(1, ev["repeats"]["conserved_bp"])
    like_all = [s for s in ev["coding"]["segments"] if s.get("coding_like")]
    coding_like = [s for s in like_all if not s.get("in_repeat")]
    transposon_orfs = [s for s in like_all if s.get("in_repeat")]
    parser_exon = ev["parser"]["exons_on_conserved"] > 0
    scored = ev["coding"]["segments"]
    para = ev.get("paralogue") or {}
    para_copy = (
        bool(like_all or parser_exon)
        and para.get("share", 0.0) >= PARALOGUE_MIN_SHARE
        and para.get("query_kmers", 0) >= PARALOGUE_MIN_KMERS
    )
    rate = ev["controls"].get("coding_like_rate")
    tail = binomial_tail(len(coding_like), len(scored), rate) if rate is not None else None
    coding_strong = tail is not None and len(coding_like) > 0 and tail <= CODING_TAIL_P
    joins = ((ev["parser"].get("best_on_conserved") or {}).get("joins_annotated_exons_of")) or []
    footing = ev["human_axis_footing"]
    borrowed = footing["kilobases_with_exon_bases"] == footing["kilobases"] and footing["kilobases"] > 0

    if borrowed:
        notes.append(
            f"every one of the {footing['kilobases']} Gnocchi kilobases it touches "
            "holds exon bases of a neighbour "
            f"({footing['exon_bp_in_kilobases']} bp): the human axis may be borrowed"
        )

    if ages.get("structured_rna", 0) / cons_bp >= 0.5:
        cls, conf = "structured_rna_copy", 0.5
        basis.append(
            f"repeats: {ages['structured_rna']} of {cons_bp} conserved bp lie in a structured-RNA copy "
            f"({', '.join(ev['repeats']['conserved_families'][:2])}); held because its paralogues are"
        )
    elif ages.get("primate_specific", 0) / cons_bp >= 0.5:
        cls, conf = "alignment_artefact", 0.4
        basis.append(
            f"repeats: {ages['primate_specific']} of {cons_bp} conserved bp lie in primate-specific repeats, "
            "which placental mammals cannot share: the constraint does not fit"
        )
    elif para_copy:
        cls, conf = "coding_copy", 0.45
        basis.append(
            f"coding: the {ev['paralogue']['query']} translates to a copy of {ev['paralogue']['protein']} "
            f"({para['shared_kmers']} of {para['query_kmers']} peptide {KMER}-mers shared): "
            "an unannotated pseudogene or retrocopy, whose constraint may be its parent's"
        )
    elif coding_strong:
        cls, conf = "coding_exon", 0.45 if parser_exon else 0.3
    elif parser_exon and joins:
        cls, conf = "coding_exon", 0.3
        basis.append(
            f"coding: the segment parser joins conserved sequence in the block to annotated exons of "
            f"{', '.join(joins)} on the same strand: "
            f"an unannotated exon of {joins[0]} is the simplest reading"
        )
    elif promoter:
        cls = "promoter"
        conf = 0.35 + (READER_WEIGHT if reader_ok else 0.0) + (MEASURED_WEIGHT if mpra_active else 0.0)
        basis.append(f"registry: {'/'.join(sorted(classes & set(PROMOTER_CLASSES)))} on the block")
    elif enhancer or insulator or reader_ok or vista_pos or mpra_active:
        cls = "regulatory"
        conf = (
            0.25
            + (REGISTRY_WEIGHT if enhancer or insulator else 0.0)
            + (READER_WEIGHT if reader_ok else 0.0)
            + (MEASURED_WEIGHT if vista_pos or mpra_active else 0.0)
            + (0.05 if ev.get("deleted_elements_on_block") else 0.0)
        )
        if enhancer or insulator:
            basis.append(
                "registry: "
                + "/".join(sorted(classes & set(ENHANCER_CLASSES + INSULATOR_CLASSES)))
                + " on the block"
            )
    elif _gene_end_segment(ev) is not None:
        cls, conf = "transcript_extension", 0.2
    else:
        cls, conf = "unexplained", 0.0

    if reader_ok:
        basis.append(
            f"reader: open on the element in {len(cells)} of {ev['reader']['cells']} cell types "
            f"(controls {expected:.2f})"
        )
    elif cells:
        notes.append(f"reader: open in {len(cells)} cell type(s), controls {expected:.2f}: not above chance")
    else:
        notes.append("reader: closed in every cell type read")
    if vista_pos:
        basis.append(f"VISTA positive: {', '.join(vista_pos)}")
    if gt["vista"] and not vista_pos:
        notes.append(f"VISTA negative: {', '.join(v['id'] for v in gt['vista'])}")
    if mpra_active:
        basis.append(f"lentiMPRA active: {len(mpra_active)} element(s)")
    elif gt["mpra"]:
        notes.append(f"lentiMPRA: {len(gt['mpra'])} element(s) tested, none active")
    if ev.get("deleted_elements_on_block"):
        named = [e["target"] for e in ev["deleted_elements_on_block"] if e["target"]]
        if named:
            basis.append(f"deletion (already scored): names {', '.join(sorted(set(named)))}")
    if gt["clinvar_pathogenic"]:
        basis.append(f"ClinVar pathogenic: {len(gt['clinvar_pathogenic'])} variant(s)")
    if gt["eqtl_genes"]:
        notes.append(f"GTEx eQTL variant(s) for {', '.join(gt['eqtl_genes'][:3])}")

    if transposon_orfs:
        notes.append(
            f"coding: {len(transposon_orfs)} coding-like segment(s) lie in transposon remnants "
            f"({', '.join(sorted({s['in_repeat'] for s in transposon_orfs}))}): "
            "old ORF fragments, not counted"
        )
    if cls == "coding_exon" and coding_strong:
        basis.append(
            f"coding: {len(coding_like)} of {len(scored)} conserved segments open-framed with codon usage "
            f"at or above the median canonical exon (binomial p {tail:.3g} at the control rate {rate:.3f})"
            + (", and the segment parser puts an exon on conserved sequence" if parser_exon else "")
            + (
                f"; no known protein shares its peptide "
                f"(best {para.get('protein')}, {para.get('share', 0):.0%})"
            )
        )
        if reader_ok or classes:
            notes.append("also open or registered: a transcribed locus is often both")
    elif coding_like:
        notes.append(
            f"coding: {len(coding_like)} segment(s) coding-like by codon usage, "
            f"parser {'agrees' if parser_exon else 'finds no exon there'}"
        )
    if scored and not any(s["open_frame"] for s in scored):
        ruled_out.append(f"coding exon: none of the {len(scored)} conserved segments has an open frame")
    if not cells and not classes:
        ruled_out.append(
            f"active element in the cells read: no DNase peak in any of the {ev['reader']['cells']} "
            "and no registry element"
        )
    if cls == "transcript_extension":
        g = _gene_end_segment(ev)
        basis.append(
            f"context: conserved sequence within {GENE_END_REACH:,} bp of {g['gene']}'s 3' end "
            f"({g['type']}), "
            "no chromatin mark: an extended 3' end is the simplest reading, unmeasured"
        )
    if cls == "unexplained":
        basis.append("no layer read here supports a class")

    if borrowed and cls not in ("alignment_artefact", "structured_rna_copy"):
        conf -= 0.1
    conf = round(max(0.0, min(CONFIDENCE_CAP, conf)), 2)
    return {
        "class": cls,
        "group": GROUP[cls],
        "confidence": conf,
        "basis": basis,
        "notes": notes,
        "ruled_out": ruled_out,
        "human_axis_borrowed": borrowed,
    }


def label(ev: dict[str, Any]) -> str:
    """One line: the class, the gene context and the node's target."""
    r = ev["reading"]
    fl = ev["flanks"]
    side = " / ".join(
        f"{g['gene']} ({g['type'].replace('_', ' ')}, {g['facing']})"
        for g in (fl.get("left"), fl.get("right"))
        if g
    )
    tgt = (ev.get("inferred_target") or {}).get("gene")
    return f"{r['class'].replace('_', ' ')} between {side}" + (f"; node target {tgt}" if tgt else "")


# ---- the set -------------------------------------------------------------------------------


def attach_paralogues(rows: list[dict[str, Any]]) -> None:
    """Coding-like segments and parser structures translated and held against the known proteome: a
    constrained open frame that is a copy of a known protein is a pseudogene or retrocopy first."""
    queries: dict[str, str] = {}
    for i, r in enumerate(rows):
        for j, s in enumerate(r["coding"]["segments"]):
            if s.get("coding_like") and len(s.get("peptide") or "") >= KMER:
                queries[f"{i}:s{j}"] = s["peptide"]
        best = r["parser"].get("best_on_conserved")
        if best and len(best.get("peptide") or "") >= KMER:
            queries[f"{i}:p"] = best["peptide"]
    found = paralogues(queries)
    for i, r in enumerate(rows):
        mine = [(k, v) for k, v in found.items() if k.split(":")[0] == str(i)]
        best = max(
            mine,
            key=lambda kv: (
                kv[1]["query_kmers"] >= PARALOGUE_MIN_KMERS,
                kv[1]["share"],
                kv[1]["shared_kmers"],
            ),
            default=None,
        )
        r["paralogue"] = (
            {"query": "parser structure" if best[0].endswith(":p") else "conserved segment", **best[1]}
            if best
            else None
        )


def binomial_tail(k: int, n: int, p: float) -> float:
    """P(X >= k) for X ~ Binomial(n, p)."""
    from math import comb

    if k <= 0:
        return 1.0
    return max(0.0, min(1.0, 1.0 - sum(comb(n, i) * p**i * (1 - p) ** (n - i) for i in range(k))))


def replicate_p(observed: float, replicates: list[float]) -> float | None:
    if not replicates:
        return None
    return round((1 + sum(1 for r in replicates if r >= observed)) / (1 + len(replicates)), 4)


def set_tests(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """The 69 against their control windows, as a set: the r-th control of every block forms replicate r."""
    n = len(rows)
    width = max((len(r["controls"]["_per_window"]) for r in rows), default=0)

    def reps(key):
        # a block with fewer control windows reuses them in turn
        return [
            sum(
                float(r["controls"]["_per_window"][k % len(r["controls"]["_per_window"])][key])
                for r in rows
                if r["controls"]["_per_window"]
            )
            / n
            for k in range(width)
        ]

    obs_cells = sum(len(r["reader"]["open_on_element"]) for r in rows) / n
    obs_ccre = sum(1 for r in rows if r["ccres"]) / n
    obs_parser = sum(1 for r in rows if r["parser"]["exons_on_conserved"]) / n
    rc, rr, rp = reps("cells_open"), reps("ccre"), reps("parser_exon_on_conserved")
    cand_scores = [s["score"] for r in rows for s in r["coding"]["segments"] if s["score"] is not None]
    cand_all = sum(r["coding"]["segments_scored"] for r in rows)
    ctrl_all = sum(r["controls"]["segments_scored"] for r in rows)
    ctrl_open_frame = sum(r["controls"]["segments_open_frame"] for r in rows)
    ctrl_like = sum(r["controls"]["segments_coding_like"] for r in rows)
    cand_like = sum(1 for r in rows for s in r["coding"]["segments"] if s.get("coding_like"))
    cells = max((r["reader"]["cells"] for r in rows), default=1)

    def pooled(side: str, num: str, den: str, per_cell: bool = False) -> float:
        get = (lambda r: r["bp"]) if side == "candidates" else (lambda r: r["controls"]["bp"])
        a = sum(get(r)[num] for r in rows)
        b = sum(get(r)[den] for r in rows) * (cells if per_cell else 1)
        return round(a / max(1, b), 4)

    return {
        "blocks": n,
        "replicates": width,
        "reader_cells_open_mean": {
            "candidates": round(obs_cells, 2),
            "controls": round(sum(rc) / max(1, len(rc)), 2),
            "p": replicate_p(obs_cells, rc),
        },
        "with_registry_element": {
            "candidates": round(obs_ccre, 3),
            "controls": round(sum(rr) / max(1, len(rr)), 3),
            "p": replicate_p(obs_ccre, rr),
        },
        "parser_exon_on_conserved": {
            "candidates": round(obs_parser, 3),
            "controls": round(sum(rp) / max(1, len(rp)), 3),
            "p": replicate_p(obs_parser, rp),
        },
        "bases": {
            "note": "pooled bases; open fractions are per cell type",
            "conserved_in_registry": {
                "candidates": pooled("candidates", "conserved_in_ccre", "conserved"),
                "controls_conserved": pooled("controls", "conserved_in_ccre", "conserved"),
                "controls_all_bases": pooled("controls", "in_ccre", "window"),
            },
            "conserved_open": {
                "candidates": pooled("candidates", "conserved_open", "conserved", True),
                "controls_conserved": pooled("controls", "conserved_open", "conserved", True),
                "controls_all_bases": pooled("controls", "open", "window", True),
            },
            "all_bases_open": {
                "candidates": pooled("candidates", "open", "window", True),
                "controls": pooled("controls", "open", "window", True),
            },
        },
        "conserved_segments_scored": {
            "candidates": cand_all,
            "candidates_open_frame": round(len(cand_scores) / max(1, cand_all), 3),
            "candidates_coding_like": cand_like,
            "controls": ctrl_all,
            "controls_open_frame": round(ctrl_open_frame / max(1, ctrl_all), 3),
            "controls_coding_like": ctrl_like,
            "controls_coding_like_fraction": round(ctrl_like / max(1, ctrl_all), 3),
            "candidates_coding_like_fraction": round(cand_like / max(1, cand_all), 3),
        },
    }


def run(
    results_dir: Path = RESULTS_DIR, progress=None, chroms: tuple[str, ...] | None = None
) -> dict[str, Any]:
    from genomeos.genome.segments import SegmentParser
    from genomeos.genome.signals import SignalSet

    t0 = time.time()
    rows = candidates(results_dir)
    if chroms:
        rows = [r for r in rows if r["chrom"] in chroms]
    signals = SignalSet.load(results_dir / SIGNALS.name)
    out: list[dict[str, Any]] = []
    sensitivity: dict[str, dict[str, Any]] = {}
    for chrom in CHROMS:
        mine = [r for r in rows if r["chrom"] == chrom]
        if not mine:
            continue
        ch = Chromosome(chrom, results_dir)
        truth = load_truth(chrom)
        parser = SegmentParser(signals, ch.codon_lo, 0.6)
        sensitivity[chrom] = ch.parser_sensitivity(parser)
        for r in sorted(mine, key=lambda x: x["start"]):
            out.append(read_block(ch, r, truth, parser))
            if progress:
                progress(f"{chrom}:{r['start']:,}-{r['end']:,} read")
        ch.close()
    attach_paralogues(out)
    scored = sum(e["controls"]["segments_scored"] for e in out)
    rate = sum(e["controls"]["segments_coding_like"] for e in out) / max(1, scored)
    for ev in out:
        ev["controls"]["coding_like_rate"] = round(rate, 4)
    for ev in out:
        ev["reading"] = reading(ev)
        ev["label"] = label(ev)
        if progress:
            progress(
                f"{ev['chrom']}:{ev['start']:,}-{ev['end']:,} "
                f"{ev['reading']['class']} {ev['reading']['confidence']}"
            )
    tests = set_tests(out)
    dist = Counter(e["reading"]["class"] for e in out)
    groups = Counter(e["reading"]["group"] for e in out)
    for e in out:
        e["controls"].pop("_per_window", None)
    return {
        "case": CASE,
        "blocks": len(out),
        "bp": sum(e["length"] for e in out),
        "by_class": {c: dist.get(c, 0) for c in CLASSES},
        "by_group": {g: groups.get(g, 0) for g in ("regulatory", "coding", "rna", "artefact", "unexplained")},
        "mean_confidence": round(sum(e["reading"]["confidence"] for e in out) / max(1, len(out)), 3),
        "set_tests": tests,
        "parser_positive_control": {
            "per_chromosome": sensitivity,
            "exons": sum(v["exons"] for v in sensitivity.values()),
            "recovered": sum(v["recovered"] for v in sensitivity.values()),
            "sensitivity": round(
                sum(v["recovered"] for v in sensitivity.values())
                / max(1, sum(v["exons"] for v in sensitivity.values())),
                3,
            ),
        },
        "candidates": out,
        "thresholds": {
            "merge_gap": MERGE_GAP,
            "coding_min_segment": CODING_MIN_SEGMENT,
            "stop_free_coverage": STOP_FREE_COVERAGE,
            "coding_like": (
                "open frame with mean codon log-odds at or above the chromosome's median "
                "canonical internal exon"
            ),
            "reader": f"open in >= {READER_MIN_CELLS} cell types and >= {READER_MIN_RATIO}x the control mean",
            "controls_per_side": CONTROLS_PER_SIDE,
            "control_step_min": CONTROL_STEP_MIN,
            "parser_flank": PARSER_FLANK,
            "gene_end_reach": GENE_END_REACH,
            "confidence_cap": CONFIDENCE_CAP,
        },
        "inputs": [
            "organised_<chrom> (candidates)",
            "gencode_v50_<chrom>.gff3.gz, <chrom>.fa (data/reference)",
            "domains_<chrom>, reader_<cell>_<chrom>, dnase_<cell>_<chrom>, ccres_<chrom>, rmsk_<chrom>",
            "phastConsElements100way_<chrom> (data/knowledge/constraint)",
            "constrained_targets_<chrom>, enhancer_targets_<chrom>",
            "enhancer_targets_all_<chrom> when complete",
            "data/knowledge: vista, mpra, clinvar, gwas, gtex",
            "signals_chr21 (segment parser PWMs)",
        ],
        "not_used": "no AlphaGenome call, no remote track: measured RNA over the blocks is the next test",
        "evidence": EVIDENCE,
        "cost": {"seconds": round(time.time() - t0, 1)},
    }


def run_and_save(results_dir: Path = RESULTS_DIR, progress=None) -> dict[str, Any]:
    out = run(results_dir, progress)
    save_result("syntax_candidates_genome_wide", out, results_dir)
    return out


def main() -> None:  # pragma: no cover - the job entry point
    out = run_and_save(progress=lambda m: print("  " + m, flush=True))
    print(json.dumps({k: out[k] for k in ("blocks", "by_class", "by_group", "set_tests")}, indent=1))
