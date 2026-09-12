# SPDX-License-Identifier: AGPL-3.0-or-later
"""Transcription-factor motifs over promoters and elements: the `requires:` lists, and the operators.

SEQUENCE-GRAMMAR.md §2 reads a promoter as a gene's dependency list: each binding site names
a factor that must be present for the block to run. This module derives those lists from
sequence with JASPAR 2026's CORE vertebrate collection (1,019 non-redundant profiles), and
then asks the question Albert put on 2026-09-12: which combinations recur across the genes
of one library, the way an operator recurs across statements. GRAMMAR-BY-COMPARISON.md
step 4.

Method. Each position frequency matrix becomes a log-odds matrix against a uniform
background (pseudocount 0.25 per cell); a hit is a window scoring at or above 85% of the
way from the matrix minimum to its maximum, on either strand. Scanning is indexed rather
than sliding: every sequence is indexed by its k-mers (k = 8, or the motif width when
shorter), each motif enumerates the k-mers of its most informative window that can still
reach the threshold, and only those positions are scored in full. A chromosome's promoters
(TSS ± 1 kb of the canonical transcript) against every profile take seconds.

Control. Every sequence is also shuffled once (every dinucleotide count kept, order lost) and scanned
the same way; a factor's enrichment is the share of real promoters with a hit over the share
of shuffled ones. A `requires:` entry is a factor with a hit at or above threshold and an
enrichment of 1.5 or more on the chromosome; without the enrichment the list would be the
alphabet, since a 6-base motif at 85% appears in most kilobases by chance.

Evidence. The matrix is `curated` (JASPAR); a hit is `predicted` (a score, not a binding
event); an operator (a combination recurring across a library) is `inferred` from counts
against expectation. None of it says a factor binds in any cell: that needs the reader.
"""

from __future__ import annotations

import math
import random
import time
import urllib.request
from array import array
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from genomeos.results import RESULTS_DIR, load_result, save_result

JASPAR_URL = "https://jaspar.elixir.no/download/data/2026/CORE/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
JASPAR_PATH = Path("data/knowledge/jaspar/core_vertebrates_2026.jaspar")
RELATIVE_SCORE = 0.85  # a hit scores at least this far from the matrix minimum to its maximum
PSEUDOCOUNT = 0.25
CORE_K = 8  # k-mer index width; motifs shorter than this use their whole width
PROMOTER_FLANK = 1000  # TSS ± this many bases
REQUIRES_ENRICHMENT = (
    1.5  # a factor enters `requires:` when real promoters carry it this much more than shuffled
)
MAX_POS = 1 << 11  # positions are packed in 11 bits: sequences are at most 2,047 bases
BASES = {"A": 0, "C": 1, "G": 2, "T": 3}
COMPLEMENT = {0: 3, 1: 2, 2: 1, 3: 0}
EVIDENCE = {
    "matrix": "curated: JASPAR 2026 CORE vertebrates, non-redundant",
    "hit": f"predicted: log-odds window at >= {RELATIVE_SCORE:.0%} of the matrix range, either strand",
    "requires": (
        "inferred: factors with a hit whose share of real promoters is "
        f">= {REQUIRES_ENRICHMENT}x the shuffled share"
    ),
    "operator": "inferred: factor combinations recurring across a library's promoters beyond expectation",
}


@dataclass
class Motif:
    id: str
    name: str
    counts: list[list[float]]  # 4 rows (A, C, G, T) × width
    pwm: list[list[float]] = field(default_factory=list)  # width × 4 log-odds
    threshold: float = 0.0
    minimum: float = 0.0
    maximum: float = 0.0
    core_offset: int = 0
    core_k: int = CORE_K

    @property
    def width(self) -> int:
        return len(self.counts[0])

    def prepare(self, relative: float = RELATIVE_SCORE) -> Motif:
        w = self.width
        self.pwm = []
        for j in range(w):
            col = [self.counts[b][j] for b in range(4)]
            n = sum(col) + 4 * PSEUDOCOUNT
            self.pwm.append([math.log2((c + PSEUDOCOUNT) / n / 0.25) for c in col])
        self.maximum = sum(max(c) for c in self.pwm)
        self.minimum = sum(min(c) for c in self.pwm)
        self.threshold = self.minimum + relative * (self.maximum - self.minimum)
        self.core_k = min(CORE_K, w)
        # the most informative window of core_k consecutive columns
        ic = [2 + sum((2 ** (c[b]) * 0.25) * c[b] for b in range(4)) for c in self.pwm]
        best, best_ic = 0, -1.0
        for o in range(w - self.core_k + 1):
            s = sum(ic[o : o + self.core_k])
            if s > best_ic:
                best, best_ic = o, s
        self.core_offset = best
        return self

    def feasible_cores(self) -> dict[int, float]:
        """Packed core k-mers whose partial score can still reach the threshold with the best of the rest."""
        rest_max = self.maximum - sum(
            max(c) for c in self.pwm[self.core_offset : self.core_offset + self.core_k]
        )
        need = self.threshold - rest_max
        cols = self.pwm[self.core_offset : self.core_offset + self.core_k]
        # branch and bound over the core, deepest columns' best remaining score as the bound
        suffix_max = [0.0] * (self.core_k + 1)
        for j in range(self.core_k - 1, -1, -1):
            suffix_max[j] = suffix_max[j + 1] + max(cols[j])
        out: dict[int, float] = {}

        def walk(j: int, code: int, score: float) -> None:
            if j == self.core_k:
                out[code] = score
                return
            for b in range(4):
                s = score + cols[j][b]
                if s + suffix_max[j + 1] >= need:
                    walk(j + 1, (code << 2) | b, s)

        walk(0, 0, 0.0)
        return out

    def score(self, codes: list[int], start: int) -> float | None:
        s = 0.0
        for j, col in enumerate(self.pwm):
            b = codes[start + j]
            if b < 0:
                return None
            s += col[b]
        return s


def parse_jaspar(text: str) -> list[Motif]:
    motifs: list[Motif] = []
    cur: Motif | None = None
    rows: list[list[float]] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if cur and len(rows) == 4:
                cur.counts = rows
                motifs.append(cur)
            parts = line[1:].split()
            cur = Motif(parts[0], parts[1] if len(parts) > 1 else parts[0], [])
            rows = []
        elif line[:1] in "ACGT" and cur is not None:
            nums = line.split("[", 1)[1].rsplit("]", 1)[0].split() if "[" in line else line.split()[1:]
            rows.append([float(x) for x in nums])
    if cur and len(rows) == 4:
        cur.counts = rows
        motifs.append(cur)
    return motifs


def load_motifs(path: Path = JASPAR_PATH, relative: float = RELATIVE_SCORE) -> list[Motif]:
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        req = urllib.request.Request(JASPAR_URL, headers={"User-Agent": "GenomeOS/0.9 (motifs)"})
        with urllib.request.urlopen(req, timeout=300) as r:  # noqa: S310
            path.write_bytes(r.read())
    return [m.prepare(relative) for m in parse_jaspar(path.read_text())]


JASPAR_TRANSFAC_URL = "https://jaspar.elixir.no/download/data/2026/CORE/JASPAR2026_CORE_vertebrates_non-redundant_pfms_transfac.txt"
JASPAR_TRANSFAC_PATH = Path("data/knowledge/jaspar/core_vertebrates_2026.transfac")
ZINC_FINGER_CLASS = "C2H2 zinc finger factors"


def load_families(path: Path = JASPAR_TRANSFAC_PATH) -> dict[str, dict[str, str]]:
    """Factor name (upper case) to its curated TFClass family and class, from JASPAR's TRANSFAC file.

    Matrices of one family share a binding mode and hit the same sites, so counting them one by one
    reads one site many times (the ZRS's 23 conserved factors are 10 families). Every record carries
    `CC tf_family:` and `CC tf_class:`; a heterodimer lists two families separated by "; ".
    """
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        req = urllib.request.Request(JASPAR_TRANSFAC_URL, headers={"User-Agent": "GenomeOS/0.9 (motifs)"})
        with urllib.request.urlopen(req, timeout=300) as r:  # noqa: S310
            path.write_bytes(r.read())
    out: dict[str, dict[str, str]] = {}
    cur: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if line.startswith("ID "):
            cur["name"] = line[3:].strip()
        elif line.startswith("CC tf_family:"):
            cur["family"] = line.split(":", 1)[1].strip()
        elif line.startswith("CC tf_class:"):
            cur["class"] = line.split(":", 1)[1].strip()
        elif line.startswith("//"):
            if cur.get("name"):
                out[cur["name"].upper()] = {"family": cur.get("family", ""), "class": cur.get("class", "")}
            cur = {}
    return out


def family_unit(factor: str, families: dict[str, dict[str, str]]) -> str:
    """The unit a factor is counted as: its TFClass family, unless the family says nothing about the motif.

    - C2H2 zinc-finger families are structural ("more than 3 adjacent zinc fingers" holds 195 matrices
      that bind unrelated sequences), so a zinc-finger factor stays itself;
    - a factor without a family annotation stays itself;
    - a heterodimer's families are joined in a fixed order, so ETV5::FOXO1 is "Ets-related + FOX".
    """
    rec = families.get(factor.upper())
    if not rec or not rec.get("family"):
        return factor.upper()
    if ZINC_FINGER_CLASS in rec.get("class", ""):
        return factor.upper()
    parts = sorted({p.strip() for p in rec["family"].split(";") if p.strip()})
    return " + ".join(parts)


def gc_content(seq: str) -> float | None:
    s = seq.upper()
    acgt = sum(s.count(b) for b in "ACGT")
    return round((s.count("G") + s.count("C")) / acgt, 4) if acgt else None


INTERSPERSED = {"SINE", "LINE", "LTR", "DNA", "Retroposon"}


def repeat_fraction(intervals: list[tuple[int, int]], repeats: list[tuple[int, int]]) -> list[float]:
    """For each interval, the share of its bases covered by the (sorted, possibly overlapping) repeats."""
    import bisect

    merged: list[tuple[int, int]] = []
    for a, b in sorted(repeats):
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a, b))
    starts = [a for a, _ in merged]
    out = []
    for lo, hi in intervals:
        k = max(0, bisect.bisect_right(starts, lo) - 1)
        cov = 0
        for a, b in merged[k:]:
            if a >= hi:
                break
            if b > lo:
                cov += min(b, hi) - max(a, lo)
        out.append(round(cov / (hi - lo), 4) if hi > lo else 0.0)
    return out


def promoter_composition(chrom: str) -> dict[str, dict[str, float]]:
    """GC and interspersed-repeat share of every canonical promoter (TSS ± flank) on a chromosome.

    Both confound a motif-pair statistic: GC-rich profiles co-hit GC-rich promoters, and profiles
    that match a transposon co-hit every promoter carrying one (promoters holding both ZNF135 and
    ZNF460 sites are 19% Alu against 1% for neither). The repeat share is read from RepeatMasker.
    """
    from genomeos.genome.repeats import load_repeats

    prom = promoters(chrom)
    loci = promoter_loci(chrom)
    if prom is None or loci is None:
        return {}
    names, seqs = prom
    where = {sym: (a, b) for sym, a, b in loci}
    reps = [(r.start, r.end) for r in load_repeats(chrom) if r.cls in INTERSPERSED]
    ivs = [where[n] for n in names]
    fracs = repeat_fraction(ivs, reps) if reps else [None] * len(names)
    out: dict[str, dict[str, float]] = {}
    for n, seq, rf in zip(names, seqs, fracs, strict=True):
        g = gc_content(seq)
        if g is not None:
            out[n] = {"gc": g, "repeat": rf}
    return out


def encode(seq: str) -> list[int]:
    return [BASES.get(c, -1) for c in seq.upper()]


def reverse_codes(codes: list[int]) -> list[int]:
    return [COMPLEMENT.get(b, -1) for b in reversed(codes)]


class Index:
    """k-mer → packed (sequence, strand, position) for a set of sequences, both strands, for several k."""

    def __init__(self, seqs: list[str], ks: set[int]):
        if any(len(s) > MAX_POS for s in seqs):
            raise ValueError(f"sequences must be shorter than {MAX_POS} bases")
        self.fwd = [encode(s) for s in seqs]
        self.rev = [reverse_codes(c) for c in self.fwd]
        self.tables: dict[int, dict[int, array]] = {}
        for k in ks:
            table: dict[int, array] = defaultdict(lambda: array("I"))
            mask = (1 << (2 * k)) - 1
            for i, (f, r) in enumerate(zip(self.fwd, self.rev, strict=True)):
                for strand, codes in ((0, f), (1, r)):
                    code, valid = 0, 0
                    for pos, b in enumerate(codes):
                        if b < 0:
                            valid = 0
                            code = 0
                            continue
                        code = ((code << 2) | b) & mask
                        valid += 1
                        if valid >= k:
                            table[code].append((i << 12) | (strand << 11) | (pos - k + 1))
            self.tables[k] = table

    def hits(self, motif: Motif) -> dict[int, tuple[float, int, int]]:
        """Best hit per sequence: sequence index → (relative score, forward start, strand)."""
        table = self.tables[motif.core_k]
        w, span = motif.width, motif.maximum - motif.minimum
        best: dict[int, tuple[float, int, int]] = {}
        for code in motif.feasible_cores():
            for packed in table.get(code, ()):
                i, strand, pos = packed >> 12, (packed >> 11) & 1, packed & (MAX_POS - 1)
                start = pos - motif.core_offset
                codes = self.rev[i] if strand else self.fwd[i]
                if start < 0 or start + w > len(codes):
                    continue
                s = motif.score(codes, start)
                if s is None or s < motif.threshold:
                    continue
                rel = (s - motif.minimum) / span if span else 1.0
                fwd_start = len(codes) - start - w if strand else start
                if i not in best or rel > best[i][0]:
                    best[i] = (round(rel, 3), fwd_start, strand)
        return best


def scan(motifs: list[Motif], seqs: list[str]) -> dict[str, dict[int, tuple[float, int, int]]]:
    """Factor name → {sequence index: best hit}; motifs of one factor (several profiles) are merged."""
    idx = Index(seqs, {m.core_k for m in motifs})
    out: dict[str, dict[int, tuple[float, int, int]]] = {}
    for m in motifs:
        h = idx.hits(m)
        if not h:
            continue
        d = out.setdefault(m.name.upper(), {})
        for i, hit in h.items():
            if i not in d or hit[0] > d[i][0]:
                d[i] = hit
    return out


def dinucleotide_shuffle(seq: str, rng: random.Random) -> str:
    """A random sequence with the same dinucleotide counts (Altschul and Erickson 1985).

    Single-base shuffling keeps only the composition and destroys the CpG and GC runs real
    promoters have, which lets long GC-rich matrices read as enriched against nothing.
    Keeping every dinucleotide count keeps that local structure in the control.
    """
    if len(seq) < 3:
        return seq
    edges: dict[str, list[str]] = defaultdict(list)
    for a, b in zip(seq, seq[1:], strict=False):
        edges[a].append(b)
    last = seq[-1]
    letters = list(edges)
    for _ in range(100):  # draw the last edge of every vertex until they form a path tree into `last`
        last_edge: dict[str, str] = {}
        for v in letters:
            if v != last:
                last_edge[v] = rng.choice(edges[v])
        ok = True
        for v in letters:
            if v == last:
                continue
            seen, u = set(), v
            while u != last and u in last_edge and u not in seen:
                seen.add(u)
                u = last_edge[u]
            if u != last:
                ok = False
                break
        if ok:
            break
    else:  # pragma: no cover - degenerate sequences
        return seq
    pools: dict[str, list[str]] = {}
    for v in letters:
        rest = list(edges[v])
        if v != last:
            rest.remove(last_edge[v])
        rng.shuffle(rest)
        if v != last:
            rest.append(last_edge[v])
        pools[v] = rest
    out = [seq[0]]
    u = seq[0]
    for _ in range(len(seq) - 1):
        nxt = pools[u].pop(0)
        out.append(nxt)
        u = nxt
    return "".join(out)


def shuffled(seqs: list[str], seed: int = 1) -> list[str]:
    rng = random.Random(seed)
    return [dinucleotide_shuffle(s, rng) for s in seqs]


def enrichment(real: dict[str, dict], control: dict[str, dict], n: int) -> dict[str, dict]:
    """Per factor: share of real sequences with a hit, share of shuffled ones, and the ratio."""
    out = {}
    for tf in set(real) | set(control):
        a, b = len(real.get(tf, {})), len(control.get(tf, {}))
        out[tf] = {
            "sequences": a,
            "shuffled": b,
            "share": round(a / n, 4) if n else None,
            "enrichment": round((a + 0.5) / (b + 0.5), 2),
        }
    return out


def promoter_loci(chrom: str, flank: int = PROMOTER_FLANK) -> list[tuple[str, int, int]] | None:
    """(symbol, start, end) of every canonical promoter, TSS ± flank; the one place the TSS is chosen."""
    from genomeos.coords import Strand
    from genomeos.genome.annotation import Annotation, default_gencode

    gff = default_gencode({chrom})
    if gff is None:
        return None
    ann = Annotation.from_gff3(gff, {chrom})
    out: list[tuple[str, int, int]] = []
    for g in ann.protein_coding():
        if g.locus.chrom != chrom:
            continue
        ts = list(g.transcripts.values())
        canon = [t for t in ts if "Ensembl_canonical" in t.tags] or sorted(
            ts, key=lambda t: -(t.locus.end - t.locus.start)
        )
        if not canon:
            continue
        t = canon[0]
        tss = t.locus.start if t.locus.strand == Strand.PLUS else t.locus.end
        out.append((g.symbol, max(0, tss - flank), tss + flank))
    return out


def promoters(chrom: str, flank: int = PROMOTER_FLANK) -> tuple[list[str], list[str]] | None:
    """Gene symbols and their promoter sequences (TSS ± flank of the canonical transcript)."""
    from genomeos.coords import Locus
    from genomeos.genome.fetch import REFERENCE
    from genomeos.genome.genome import Genome

    loci = promoter_loci(chrom, flank)
    fa = REFERENCE / f"{chrom}.fa.gz"
    if loci is None or not fa.exists():
        return None
    genome = Genome.from_fasta(fa)
    names, seqs = [], []
    for sym, start, end in loci:
        seq = str(genome.fetch(Locus(chrom, start, end)))
        if len(seq) < 2 * flank:
            continue
        names.append(sym)
        seqs.append(seq)
    return names, seqs


def elements(chrom: str, results_dir: Path = RESULTS_DIR) -> tuple[list[dict], list[str]] | None:
    from genomeos.coords import Locus
    from genomeos.genome.fetch import REFERENCE
    from genomeos.genome.genome import Genome

    fa = REFERENCE / f"{chrom}.fa.gz"
    if not fa.exists():
        return None
    seen: set[str] = set()
    rows: list[dict] = []
    for name in ("constrained_targets", "enhancer_targets"):
        r = load_result(f"{name}_{chrom}", results_dir) or {}
        for e in r.get("elements", []):
            if e["id"] in seen or e["end"] - e["start"] >= MAX_POS:
                continue
            seen.add(e["id"])
            pc = e.get("predicted_coding") or {}
            rows.append({"id": e["id"], "start": e["start"], "end": e["end"], "target": pc.get("gene")})
    genome = Genome.from_fasta(fa)
    seqs = [str(genome.fetch(Locus(chrom, e["start"], e["end"]))) for e in rows]
    return rows, seqs


def requires(
    hits_for_seq: dict[str, tuple[float, int, int]], enrich: dict[str, dict], top: int = 8
) -> list[dict]:
    """A sequence's `requires:` list: enriched factors with a hit, best score first."""
    rows = [
        {
            "factor": tf,
            "score": hit[0],
            "position": hit[1],
            "strand": "-" if hit[2] else "+",
            "enrichment": enrich[tf]["enrichment"],
        }
        for tf, hit in hits_for_seq.items()
        if enrich.get(tf, {}).get("enrichment", 0) >= REQUIRES_ENRICHMENT
    ]
    rows.sort(key=lambda r: (-r["score"], -r["enrichment"]))
    return rows[:top]


def library_members() -> dict[str, list[str]] | None:
    try:
        from genomeos.lib.membership import KnowledgeBase

        d = KnowledgeBase.distilled()
        return d.get("members") if d else None
    except Exception:  # pragma: no cover - optional
        return None


def library_factors(
    gene_hits: dict[str, dict[str, float]], members: dict[str, list[str]], min_members: int = 5
) -> dict[str, dict]:
    """Per library: factors hitting a larger share of its members' promoters than of all promoters."""
    n = len(gene_hits)
    if not n:
        return {}
    all_share = Counter()
    for hits in gene_hits.values():
        all_share.update(hits.keys())
    out: dict[str, dict] = {}
    for lib, syms in members.items():
        here = [gene_hits[s] for s in syms if s in gene_hits]
        if len(here) < min_members:
            continue
        cnt = Counter()
        for h in here:
            cnt.update(h.keys())
        rows = []
        for tf, c in cnt.items():
            share, base = c / len(here), all_share[tf] / n
            if c >= 3 and base > 0 and share / base >= 1.5:
                rows.append(
                    {
                        "factor": tf,
                        "members_with_hit": c,
                        "share": round(share, 3),
                        "all_share": round(base, 3),
                        "ratio": round(share / base, 2),
                    }
                )
        rows.sort(key=lambda r: (-r["ratio"], -r["members_with_hit"]))
        out[lib] = {"members_on_chromosome": len(here), "factors": rows[:12]}
    return out


def build(
    chrom: str, results_dir: Path = RESULTS_DIR, motifs: list[Motif] | None = None, progress=None
) -> dict:
    t0 = time.time()
    motifs = motifs or load_motifs()
    prom = promoters(chrom)
    if prom is None:
        raise FileNotFoundError(
            f"no sequence or GENCODE models for {chrom}; run genomeos data fetch --chrom {chrom}"
        )
    names, seqs = prom
    real = scan(motifs, seqs)
    ctrl = scan(motifs, shuffled(seqs))
    enrich = enrichment(real, ctrl, len(seqs))
    if progress:
        progress(
            f"{chrom}: {len(seqs)} promoters scanned against {len(motifs)} profiles "
            f"in {time.time() - t0:.0f} s"
        )
    genes: dict[str, dict] = {}
    gene_hits: dict[str, dict[str, float]] = {}
    for i, sym in enumerate(names):
        hits_i = {tf: h[i] for tf, h in real.items() if i in h}
        req = requires(hits_i, enrich)
        genes[sym] = {"factors_hit": len(hits_i), "requires": req}
        gene_hits[sym] = {r["factor"]: r["score"] for r in req}
    el = elements(chrom, results_dir)
    element_rows: list[dict] = []
    if el and el[1]:
        rows, eseqs = el
        ereal = scan(motifs, eseqs)
        for i, e in enumerate(rows):
            hits_i = {tf: h[i] for tf, h in ereal.items() if i in h}
            element_rows.append({**e, "factors_hit": len(hits_i), "requires": requires(hits_i, enrich)})
    members = library_members() or {}
    ranked = sorted(enrich.items(), key=lambda kv: (-kv[1]["enrichment"], -kv[1]["sequences"]))
    return {
        "chrom": chrom,
        "profiles": len(motifs),
        "factors": len({m.name.upper() for m in motifs}),
        "promoters": len(seqs),
        "flank": PROMOTER_FLANK,
        "thresholds": {"relative_score": RELATIVE_SCORE, "requires_enrichment": REQUIRES_ENRICHMENT},
        "genes": genes,
        "elements": element_rows,
        "factor_enrichment": {tf: v for tf, v in ranked},
        "most_enriched": [{"factor": tf, **v} for tf, v in ranked[:25] if v["sequences"] >= 3],
        "requires_per_promoter": round(sum(len(g["requires"]) for g in genes.values()) / len(genes), 2)
        if genes
        else None,
        "libraries": library_factors(gene_hits, members),
        "evidence": EVIDENCE,
        "cost": {"seconds": round(time.time() - t0, 1)},
    }


def run_and_save(
    chrom: str, results_dir: Path = RESULTS_DIR, motifs: list[Motif] | None = None, progress=None
) -> dict:
    out = build(chrom, results_dir, motifs, progress)
    save_result(f"motifs_{chrom}", out, results_dir)
    return out


OPERATOR_RATIO = 2.0  # a pair's observed co-occurrence over its GC-matched expectation
OPERATOR_SHARE = 0.25  # and present together in at least this share of the library's members
GC_BINS = 10
REPEAT_EDGES = (0.0, 0.1, 0.3)  # repeat share bins: none, up to 10%, up to 30%, more
MIN_STRATUM = 30  # a GC-by-repeat stratum thinner than this falls back to its GC decile


def repeat_bin(share: float | None) -> int | None:
    if share is None:
        return None
    if share <= REPEAT_EDGES[0]:
        return 0
    return 1 if share <= REPEAT_EDGES[1] else 2 if share <= REPEAT_EDGES[2] else 3


def gc_bins(gc: dict[str, float], bins: int = GC_BINS) -> dict[str, int]:
    """Gene to GC decile (0 lowest); genes with no GC value are left out."""
    ranked = sorted(gc.items(), key=lambda kv: kv[1])
    n = len(ranked)
    return {g: min(bins - 1, i * bins // n) for i, (g, _) in enumerate(ranked)} if n else {}


def distil(
    results_dir: Path = RESULTS_DIR,
    members: dict[str, list[str]] | None = None,
    min_members: int = 8,
    families: dict[str, dict[str, str]] | None = None,
    composition: dict[str, dict[str, float]] | None = None,
) -> dict:
    """The genome: `requires:` per gene pooled, and per library the enriched factor families and the family
    pairs that recur across its members' promoters beyond a GC-matched expectation (operators).

    Two corrections over counting matrices against global shares. Units are TFClass families
    (`family_unit`), so near-identical matrices (MEF2A and MEF2D) are one unit and cannot form a pair.
    Expectations are taken within strata of GC decile and interspersed-repeat share: a member is
    expected to carry the units promoters of its own composition carry, so GC-rich profiles co-hitting
    (CGGBP1 with ZNF93) and transposon-matching profiles co-hitting (ZNF135 with ZNF460, on Alu) are
    explained by composition instead of reported as logic. The naive expectation is kept beside it.
    """
    import glob
    from itertools import combinations

    families = families if families is not None else load_families()
    if composition is None:
        composition = (load_result("promoter_composition_genome_wide", results_dir) or {}).get("genes") or {}
    gene_units: dict[str, set[str]] = {}
    chroms = 0
    seconds = 0.0
    for f in sorted(glob.glob(str(results_dir / "motifs_chr*.json"))):
        r = load_result(Path(f).stem, results_dir)
        if not r or "genes" not in r:
            continue
        chroms += 1
        seconds += (r.get("cost") or {}).get("seconds", 0)
        for sym, g in r["genes"].items():
            gene_units[sym] = {family_unit(x["factor"], families) for x in g.get("requires", [])}
    n = len(gene_units)
    share = Counter()
    for us in gene_units.values():
        share.update(us)
    unit_share = {u: round(c / n, 4) for u, c in share.items()} if n else {}
    binned = gc_bins(
        {g: v["gc"] for g, v in composition.items() if g in gene_units and v.get("gc") is not None}
    )
    stratum = {g: (b, repeat_bin(composition[g].get("repeat"))) for g, b in binned.items()}
    gc_size, gc_counts = Counter(binned.values()), defaultdict(Counter)
    st_size, st_counts = Counter(stratum.values()), defaultdict(Counter)
    for g, b in binned.items():
        gc_counts[b].update(gene_units[g])
        st_counts[stratum[g]].update(gene_units[g])

    def p_unit(u: str, gene: str) -> float:
        st = stratum.get(gene)
        if st is not None and st[1] is not None and st_size[st] >= MIN_STRATUM:
            return st_counts[st][u] / st_size[st]
        b = binned.get(gene)
        if b is not None and gc_size[b]:
            return gc_counts[b][u] / gc_size[b]
        return unit_share.get(u, 0.0)

    members = members if members is not None else (library_members() or {})
    libraries: dict[str, dict] = {}
    explained_total = 0
    for lib, syms in members.items():
        here = {s_: gene_units[s_] for s_ in syms if s_ in gene_units}
        m = len(here)
        if m < min_members:
            continue
        cnt = Counter()
        for us in here.values():
            cnt.update(us)
        factors = []
        for u, c in cnt.items():
            expected = sum(p_unit(u, g) for g in here) / m
            if c >= 3 and expected and (c / m) / expected >= 1.5:
                factors.append(
                    {
                        "family": u,
                        "members_with_hit": c,
                        "share": round(c / m, 3),
                        "expected_share": round(expected, 4),
                        "all_share": unit_share.get(u, 0),
                        "ratio": round((c / m) / expected, 2),
                    }
                )
        factors.sort(key=lambda r: (-r["ratio"], -r["members_with_hit"]))
        pairs = Counter()
        for us in here.values():
            for a, b in combinations(sorted(us), 2):
                pairs[(a, b)] += 1
        ops = []
        explained = 0
        for (a, b), c in pairs.items():
            if c < 3 or c / m < OPERATOR_SHARE:
                continue
            naive = m * unit_share.get(a, 0) * unit_share.get(b, 0)
            expected = sum(p_unit(a, g) * p_unit(b, g) for g in here)
            naive_ratio = c / naive if naive else None
            ratio = c / expected if expected else None
            if ratio is not None and ratio >= OPERATOR_RATIO:
                ops.append(
                    {
                        "families": [a, b],
                        "members_with_both": c,
                        "share": round(c / m, 3),
                        "expected": round(expected, 2),
                        "ratio": round(ratio, 2),
                        "naive_ratio": round(naive_ratio, 2) if naive_ratio else None,
                    }
                )
            elif naive_ratio is not None and naive_ratio >= OPERATOR_RATIO:
                explained += 1
        explained_total += explained
        ops.sort(key=lambda r: (-r["ratio"], -r["members_with_both"]))
        libraries[lib] = {
            "members_placed": m,
            "members": len(syms),
            "families": factors[:12],
            "operators": ops[:12],
            "pairs_explained_by_gc": explained,
        }
    return {
        "chromosomes": chroms,
        "promoters": n,
        "units": len(unit_share),
        "promoters_with_composition": len(binned),
        "strata_used": sum(1 for st, k in st_size.items() if st[1] is not None and k >= MIN_STRATUM),
        "unit_share": unit_share,
        "libraries": libraries,
        "pairs_explained_by_gc": explained_total,
        "libraries_with_operators": sum(1 for v in libraries.values() if v["operators"]),
        "thresholds": {
            "relative_score": RELATIVE_SCORE,
            "requires_enrichment": REQUIRES_ENRICHMENT,
            "operator_ratio": OPERATOR_RATIO,
            "operator_share": OPERATOR_SHARE,
            "gc_bins": GC_BINS,
            "repeat_edges": list(REPEAT_EDGES),
            "min_stratum": MIN_STRATUM,
        },
        "evidence": {
            **EVIDENCE,
            "units": "curated: TFClass family per JASPAR matrix; C2H2 zinc-finger factors counted one by one",
            "null": "inferred: expectation within GC-decile by repeat-share strata; naive one kept too",
        },
        "cost": {"seconds": round(seconds, 1)},
    }
