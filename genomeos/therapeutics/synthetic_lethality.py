# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""What a loss makes indispensable: dependency tested against loss across cell lines.

The rest of this module's neighbours ask what the tumour *displays*. This one
asks the opposite question, the one a deep deletion actually poses. A
homozygously deleted gene is never a target — `docs/THERAPEUTICS.md` says so and
`genomeos/cancer/alterations.py` enforces it — and what a deletion points at is
the dependency the loss creates. That dependency is measurable, publicly, in
the same cells:

    cell lines that lost A   ->  do they need B more than the lines that kept A?

The test is a rank test, not a fold change: gene effect is not normally
distributed, the lost group is usually small, and a single hypersensitive line
would carry a mean. Every pair in the pre-declared candidate space is tested,
the one-sided p-values go through Benjamini-Hochberg, and a pair is called
recovered only at q <= 0.10 with the dependency stronger in the lines that lost
the gene.

**The candidate space and the controls are declared in this file before any
result was looked at.** `PREREGISTERED` names six pairs whose answer is known
from the literature and one of them is expected to *fail*: BRCA1/BRCA2 loss
against PARP1 dependency is the textbook synthetic-lethal pair and it is
famously weak in CRISPR knockout screens, because the clinical effect comes
from trapping PARP1 on DNA with an inhibitor rather than from deleting it. A
benchmark that reported PARP1 as recovered would be measuring its own wishes.
`shuffled_null` re-runs the whole sweep with the loss labels permuted, which is
the control for the sweep itself rather than for any one pair.

Nothing here is a therapy. A dependency is a hypothesis about a target; no
molecule, no construct and no protocol follows from it in this project.
"""

from __future__ import annotations

import json
import math
import random
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

from . import depmap
from .evidence import Evidence, prediction
from .providers import DiskCache

# --- pre-registration -------------------------------------------------------------------

#: False-discovery rate for the sweep, Benjamini-Hochberg, one-sided.
FDR = 0.10
#: A pair is not tested below these group sizes; an underpowered test is not a negative.
MIN_LOST = 5
MIN_INTACT = 20
#: Loss calls. A damaging mutation is `LikelyLoF` in DepMap's own matrix; a
#: homozygous deletion is PureCN absolute copy number below half a copy; MSI-high
#: is MSIsensor2 at or above 20, the threshold used in the WRN literature.
DAMAGING_MIN = 1
DELETION_MAX_COPIES = 0.5
MSI_HIGH = 20.0
#: Paralogue pairs are taken from the paralogue cache, closest first.
MIN_PARALOGUE_IDENTITY = 20.0
MAX_PARALOGUES_PER_GENE = 15
#: The shuffled-label null: how many permutations, and the seed.
NULL_PERMUTATIONS = 20
SEED = 20260916

#: Verdict rules, fixed before the run.
VERDICTS = {
    "recovered": f"q <= {FDR} and the dependency is stronger in the lines that lost the gene",
    "nominal": "p <= 0.05 but q above the false-discovery rate",
    "not_recovered": "neither, with enough lines to have seen an effect",
    "untestable": f"fewer than {MIN_LOST} lost or {MIN_INTACT} intact lines, or the gene is absent",
}

#: The controls. Written down before the first result was read: five pairs that
#: should come back, one that should not, and why in each case.
PREREGISTERED: tuple[dict[str, Any], ...] = (
    {
        "lost": "BRCA1",
        "loss": "damaging_mutation",
        "dependency": "PARP1",
        "expect": "weak",
        "why": (
            "the textbook synthetic-lethal pair clinically, and known to be weak in CRISPR "
            "knockout screens: PARP inhibitors trap PARP1 on DNA, which deleting PARP1 does not "
            "reproduce. A sweep that recovers this pair strongly is suspect, not impressive."
        ),
    },
    {
        "lost": "BRCA2",
        "loss": "damaging_mutation",
        "dependency": "PARP1",
        "expect": "weak",
        "why": "as BRCA1; the same mechanistic mismatch between knockout and inhibition",
    },
    {
        "lost": "ARID1A",
        "loss": "damaging_mutation",
        "dependency": "ARID1B",
        "expect": "recovered",
        "why": "the paralogue that carries the residual BAF complex when ARID1A is lost",
    },
    {
        "lost": "SMARCA4",
        "loss": "damaging_mutation",
        "dependency": "SMARCA2",
        "expect": "recovered",
        "why": "the other BAF ATPase; the basis of the SMARCA2 degrader programmes",
    },
    {
        "lost": "MTAP",
        "loss": "deletion",
        "dependency": "PRMT5",
        "expect": "recovered",
        "why": (
            "MTAP deletion, passenger to CDKN2A loss, accumulates MTA, which inhibits PRMT5 "
            "partially and leaves the cell dependent on what activity remains"
        ),
    },
    {
        "lost": "MTAP",
        "loss": "deletion",
        "dependency": "MAT2A",
        "expect": "recovered",
        "why": "the same axis one step back: SAM supply to the partially inhibited PRMT5",
    },
    {
        "lost": "MSI",
        "loss": "msi_high",
        "dependency": "WRN",
        "expect": "recovered",
        "why": (
            "expanded TA dinucleotide repeats in MSI-high cells form non-B DNA that needs the "
            "WRN helicase; the cleanest dependency in the public screens"
        ),
    },
)

CLAIM = (
    "Cell lines that lost gene A depend more on gene B than lines that kept A, for pairs drawn "
    "from a pre-declared DNA-repair and replication-stress panel plus paralogue pairs. Recovered "
    "at Benjamini-Hochberg q <= 0.10 one-sided, with the five positive controls recovered, the "
    "two PARP1 controls expected to be weak, and a permuted-label sweep finding nothing."
)

# --- the pre-declared candidate space ---------------------------------------------------

#: The panel: DNA repair and replication stress, by pathway, plus the chromatin
#: and methionine-salvage genes the paralogue controls live in. Membership is a
#: statement of what was asked, and is fixed before the sweep.
PATHWAYS: dict[str, tuple[str, ...]] = {
    "homologous_recombination": (
        "BRCA1",
        "BRCA2",
        "PALB2",
        "RAD51",
        "RAD51B",
        "RAD51C",
        "RAD51D",
        "RAD51AP1",
        "RAD52",
        "RAD54B",
        "RAD54L",
        "XRCC2",
        "XRCC3",
        "BARD1",
        "BRIP1",
        "RBBP8",
        "MRE11",
        "RAD50",
        "NBN",
        "ATM",
        "TP53BP1",
        "RIF1",
        "SHLD1",
        "SHLD2",
        "SHLD3",
        "BLM",
        "RMI1",
        "RMI2",
        "GEN1",
        "MUS81",
        "EME1",
        "SLX4",
        "HELQ",
        "SWSAP1",
        "ZSWIM7",
    ),
    "nonhomologous_end_joining": (
        "PRKDC",
        "XRCC4",
        "XRCC5",
        "XRCC6",
        "LIG4",
        "NHEJ1",
        "DCLRE1C",
        "POLL",
        "POLM",
    ),
    "single_strand_break_and_alternative_end_joining": (
        "POLQ",
        "PARP1",
        "PARP2",
        "PARP3",
        "PARG",
        "LIG1",
        "LIG3",
        "XRCC1",
        "POLB",
        "PNKP",
        "APTX",
        "APLF",
    ),
    "mismatch_repair": ("MLH1", "MLH3", "MSH2", "MSH3", "MSH6", "PMS1", "PMS2", "EXO1"),
    "base_excision_repair": (
        "OGG1",
        "MUTYH",
        "UNG",
        "SMUG1",
        "MBD4",
        "TDG",
        "NEIL1",
        "NEIL2",
        "NEIL3",
        "NTHL1",
        "MPG",
        "APEX1",
        "APEX2",
    ),
    "nucleotide_excision_repair": (
        "ERCC1",
        "ERCC2",
        "ERCC3",
        "ERCC4",
        "ERCC5",
        "ERCC6",
        "ERCC8",
        "XPA",
        "XPC",
        "DDB1",
        "DDB2",
        "GTF2H1",
        "GTF2H5",
        "RAD23A",
        "RAD23B",
        "CETN2",
    ),
    "fanconi_anaemia": (
        "FANCA",
        "FANCB",
        "FANCC",
        "FANCD2",
        "FANCE",
        "FANCF",
        "FANCG",
        "FANCI",
        "FANCL",
        "FANCM",
        "UBE2T",
        "RFWD3",
        "CENPS",
        "CENPX",
    ),
    "translesion_synthesis": (
        "POLH",
        "POLI",
        "POLK",
        "POLN",
        "REV1",
        "REV3L",
        "MAD2L2",
        "RAD18",
        "UBE2A",
        "UBE2B",
        "PCNA",
    ),
    "checkpoint_and_replication_stress": (
        "ATR",
        "ATRIP",
        "CHEK1",
        "CHEK2",
        "CLSPN",
        "TOPBP1",
        "RAD17",
        "RAD9A",
        "RAD1",
        "HUS1",
        "WEE1",
        "CDC25A",
        "CDC7",
        "DBF4",
        "TIMELESS",
        "TIPIN",
        "RPA1",
        "RPA2",
        "RPA3",
        "SMARCAL1",
        "ZRANB3",
        "HLTF",
        "WRNIP1",
        "ATAD5",
        "USP1",
        "WDR48",
        "SPRTN",
        "FEN1",
        "DNA2",
        "RECQL",
        "RECQL4",
        "RECQL5",
        "WRN",
        "PIF1",
        "SLFN11",
        "RNASEH2A",
        "RNASEH2B",
        "RNASEH2C",
        "SAMHD1",
        "SETX",
        "TDP1",
        "TDP2",
    ),
    "topoisomerases": ("TOP1", "TOP2A", "TOP2B", "TOP3A", "TOP3B", "ZNF451"),
    "chromatin_paralogue_pairs": (
        "ARID1A",
        "ARID1B",
        "ARID2",
        "SMARCA2",
        "SMARCA4",
        "SMARCB1",
        "PBRM1",
        "SMARCC1",
        "SMARCC2",
        "SMARCE1",
        "BRD7",
        "BRD9",
        "EZH2",
        "EED",
        "SUZ12",
        "KDM6A",
        "UTY",
        "KDM6B",
        "CREBBP",
        "EP300",
        "STAG1",
        "STAG2",
        "VPS4A",
        "VPS4B",
        "ASF1A",
        "ASF1B",
        "DNMT1",
        "DNMT3A",
        "DNMT3B",
        "KMT2C",
        "KMT2D",
        "KAT6A",
        "KAT6B",
    ),
    "methionine_salvage": ("MTAP", "PRMT5", "MAT2A", "WDR77", "RIOK1", "CDKN2A", "CDKN2B"),
}

PANEL: tuple[str, ...] = tuple(sorted({g for genes in PATHWAYS.values() for g in genes}))

ENSEMBL_SYMBOL_LOOKUP = "https://rest.ensembl.org/lookup/symbol/homo_sapiens"

# --- statistics -------------------------------------------------------------------------


def _phi(z: float) -> float:
    """Standard normal CDF."""
    return 0.5 * math.erfc(-z / math.sqrt(2.0))


def ranks(values: dict[str, float]) -> tuple[dict[str, float], float]:
    """Average ranks of a model -> value map, with the tie correction term.

    Returned once per (dependency gene, coverage set) and reused for every loss
    gene sharing that coverage, which is what makes a permutation null cheap.
    """
    ordered = sorted(values.items(), key=lambda kv: kv[1])
    out: dict[str, float] = {}
    ties = 0.0
    i = 0
    while i < len(ordered):
        j = i
        while j + 1 < len(ordered) and ordered[j + 1][1] == ordered[i][1]:
            j += 1
        rank = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            out[ordered[k][0]] = rank
        run = j - i + 1
        if run > 1:
            ties += run**3 - run
        i = j + 1
    return out, ties


def rank_sum_p(rank_of: dict[str, float], ties: float, group: list[str]) -> tuple[float, float, float]:
    """One-sided Mann-Whitney: is `group` stochastically *lower* than the rest?

    Returns (p, Cliff's delta, U). Lower gene effect is a stronger dependency,
    so a negative delta is the direction synthetic lethality predicts. The
    normal approximation carries the tie correction and a continuity correction;
    with the group sizes this sweep requires it is the right approximation and
    is stated as such rather than called exact.
    """
    n = len(rank_of)
    n_a = len(group)
    n_b = n - n_a
    if n_a == 0 or n_b == 0:
        return 1.0, 0.0, 0.0
    r_a = math.fsum(rank_of[m] for m in group)
    u_a = r_a - n_a * (n_a + 1) / 2.0
    mean = n_a * n_b / 2.0
    var = (n_a * n_b / (n * (n - 1.0))) * ((n**3 - n) / 12.0 - ties / 12.0)
    delta = 2.0 * u_a / (n_a * n_b) - 1.0
    if var <= 0:
        return 1.0, delta, u_a
    z = (u_a - mean + 0.5) / math.sqrt(var)
    return _phi(z), delta, u_a


def benjamini_hochberg(pvalues: list[float]) -> list[float]:
    """BH q-values, in the order the p-values came in."""
    n = len(pvalues)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: pvalues[i])
    q = [1.0] * n
    running = 1.0
    for rank, i in reversed(list(enumerate(order, start=1))):
        running = min(running, pvalues[i] * n / rank)
        q[i] = min(1.0, running)
    return q


def median(values: list[float]) -> float | None:
    if not values:
        return None
    s = sorted(values)
    mid = len(s) // 2
    return s[mid] if len(s) % 2 else (s[mid - 1] + s[mid]) / 2.0


# --- loss calls -------------------------------------------------------------------------


def loss_status(
    gene: str,
    kind: str,
    damaging: dict[str, dict[str, float]],
    absolute_cn: dict[str, dict[str, float]],
    msi: dict[str, float],
) -> dict[str, bool]:
    """Model -> lost (True) or kept (False); models with no call are absent.

    Four definitions, each named in the result rather than blended:

    * `damaging_mutation` - DepMap's own likely-loss-of-function call
    * `deletion` - PureCN absolute copy number below half a copy
    * `loss` - either of those; a model counts as *kept* only when both calls
      are present and neither fires, because a missing copy-number value cannot
      rule out a deletion
    * `msi_high` - not a gene at all but the MSIsensor2 state of the line
    """
    if kind == "msi_high":
        return {m: v >= MSI_HIGH for m, v in msi.items()}
    out: dict[str, bool] = {}
    if kind in ("damaging_mutation", "loss"):
        for m, row in damaging.items():
            if gene in row:
                out[m] = row[gene] >= DAMAGING_MIN
    if kind == "deletion":
        return {m: row[gene] < DELETION_MAX_COPIES for m, row in absolute_cn.items() if gene in row}
    if kind == "loss":
        merged: dict[str, bool] = {}
        for m, lost in out.items():
            cn = absolute_cn.get(m, {}).get(gene)
            if lost:
                merged[m] = True
            elif cn is not None:
                merged[m] = cn < DELETION_MAX_COPIES
        return merged
    return out


# --- the sweep --------------------------------------------------------------------------


class _RankStore:
    """Ranks per (dependency gene, coverage set), computed once and reused."""

    def __init__(self, effect: dict[str, dict[str, float]]) -> None:
        self.effect = effect
        self._cache: dict[tuple[str, int], tuple[dict[str, float], float]] = {}
        self._sets: dict[frozenset[str], frozenset[str]] = {}

    def canonical(self, models: frozenset[str]) -> frozenset[str]:
        return self._sets.setdefault(models, models)

    def get(self, gene: str, coverage: frozenset[str]) -> tuple[dict[str, float], float] | None:
        key = (gene, id(coverage))
        hit = self._cache.get(key)
        if hit is not None:
            return hit
        values = {m: row[gene] for m, row in self.effect.items() if gene in row and m in coverage}
        if not values:
            return None
        got = ranks(values)
        self._cache[key] = got
        return got


def test_pair(
    lost_gene: str,
    dependency: str,
    kind: str,
    status: dict[str, bool],
    store: _RankStore,
    coverage: frozenset[str] | None = None,
    min_lost: int = MIN_LOST,
    min_intact: int = MIN_INTACT,
) -> dict[str, Any]:
    """One pair: do the lines that lost A depend more on B?"""
    coverage = coverage if coverage is not None else store.canonical(frozenset(status))
    got = store.get(dependency, coverage)
    row: dict[str, Any] = {
        "lost": lost_gene,
        "loss": kind,
        "dependency": dependency,
    }
    if got is None:
        return {**row, "verdict": "untestable", "reason": f"no gene effect for {dependency}"}
    rank_of, ties = got
    lost = [m for m in rank_of if status.get(m)]
    intact = [m for m in rank_of if status.get(m) is False]
    row.update({"lines_lost": len(lost), "lines_intact": len(intact)})
    if len(lost) < min_lost or len(intact) < min_intact:
        return {
            **row,
            "verdict": "untestable",
            "reason": f"{len(lost)} lost and {len(intact)} intact lines with a {dependency} gene effect",
        }
    p, delta, _u = rank_sum_p(rank_of, ties, lost)
    effect = store.effect
    lost_values = [effect[m][dependency] for m in lost]
    intact_values = [effect[m][dependency] for m in intact]
    m_lost = median(lost_values)
    m_intact = median(intact_values)
    return {
        **row,
        "median_effect_lost": round(m_lost, 3) if m_lost is not None else None,
        "median_effect_intact": round(m_intact, 3) if m_intact is not None else None,
        "median_difference": round((m_lost or 0.0) - (m_intact or 0.0), 3),
        "cliffs_delta": round(delta, 3),
        "p": p,
    }


def candidate_pairs(
    paralogues: dict[str, list[dict[str, Any]]],
    panel: tuple[str, ...] = PANEL,
) -> tuple[list[tuple[str, str]], dict[str, Any]]:
    """The pre-declared candidate space: panel x panel, plus paralogue pairs.

    Paralogue pairs come from the project's own paralogue cache
    (`data/knowledge/therapeutics/paralogues`, Ensembl Compara), both
    directions, closest first, above a declared identity floor.
    """
    pairs: set[tuple[str, str]] = set()
    for a in panel:
        for b in panel:
            if a != b:
                pairs.add((a, b))
    para_pairs: set[tuple[str, str]] = set()
    for gene, rows in paralogues.items():
        kept = [
            r
            for r in rows
            if (r.get("identity") or 0) >= MIN_PARALOGUE_IDENTITY and (r.get("symbol") or "") != gene
        ][:MAX_PARALOGUES_PER_GENE]
        for r in kept:
            partner = (r.get("symbol") or "").upper()
            if not partner or partner == gene:
                continue
            para_pairs.add((gene, partner))
            para_pairs.add((partner, gene))
    for control in PREREGISTERED:
        if control["lost"] != "MSI":
            pairs.add((control["lost"], control["dependency"]))
    space = {
        "panel_genes": len(panel),
        "pathways": {k: len(v) for k, v in PATHWAYS.items()},
        "panel_pairs": len(pairs),
        "paralogue_pairs": len(para_pairs - pairs),
        "paralogue_identity_floor": MIN_PARALOGUE_IDENTITY,
        "paralogues_per_gene_cap": MAX_PARALOGUES_PER_GENE,
    }
    return sorted(pairs | para_pairs), space


# --- paralogues, through the project's existing cache ------------------------------------


def ensembl_ids(genes: list[str], net: bool = True, log=None) -> dict[str, str]:
    """Symbol -> Ensembl gene id, cached; the paralogue lookup is keyed by id."""
    cache = DiskCache("ensembl_gene_ids")
    out: dict[str, str] = {}
    missing: list[str] = []
    for g in genes:
        hit = cache.get(g)
        if hit:
            out[g] = hit
        else:
            missing.append(g)
    if missing and net:
        for i in range(0, len(missing), 200):
            chunk = missing[i : i + 200]
            body = json.dumps({"symbols": chunk}).encode()
            req = urllib.request.Request(
                f"{ENSEMBL_SYMBOL_LOOKUP}?content-type=application/json",
                data=body,
                headers={
                    "Content-Type": "application/json",
                    "Accept": "application/json",
                    "User-Agent": "GenomeOS/0.1 (therapeutics)",
                },
            )
            try:
                with urllib.request.urlopen(req, timeout=60) as r:  # noqa: S310
                    res = json.load(r)
            except (urllib.error.URLError, OSError, ValueError) as e:
                if log:
                    print(f"  Ensembl symbol lookup failed: {type(e).__name__}: {e}", file=log)
                break
            for symbol, rec in (res or {}).items():
                ident = (rec or {}).get("id")
                if ident:
                    out[symbol] = ident
                    cache.put(symbol, ident)
    return out


def paralogue_table(genes: list[str], net: bool = True, log=None) -> dict[str, list[dict[str, Any]]]:
    """Paralogues of each gene, from the project's paralogue cache, filled from Ensembl."""
    from .providers import EnsemblParalogueProvider

    provider = EnsemblParalogueProvider(net=net)
    ids = ensembl_ids(genes, net=net, log=log)
    out: dict[str, list[dict[str, Any]]] = {}
    for gene in genes:
        answer = provider.paralogues(gene, ids.get(gene))
        if not answer.available:
            continue
        out[gene] = [
            {"symbol": (r.get("symbol") or "").upper(), "identity": r.get("identity")}
            for r in (answer.data or [])
        ]
    if log:
        print(f"  paralogues: {len(out)}/{len(genes)} genes resolved", file=log, flush=True)
    return out


# --- the run ----------------------------------------------------------------------------


#: How much two loss calls have to overlap before the sweep reports them as one event.
SAME_EVENT_OVERLAP = 0.8


def co_loss(
    statuses: dict[str, dict[str, bool]],
    genes: list[str],
    threshold: float = SAME_EVENT_OVERLAP,
) -> list[dict[str, Any]]:
    """Loss genes lost in the same lines: a pair test cannot tell them apart.

    A co-deleted neighbour is the commonest way a dependency test is right about
    the statistics and wrong about the mechanism - MTAP is deleted with CDKN2A at
    9p21, so a dependency on the methylosome appears under both names. The sweep
    cannot resolve which gene carries the effect, and this is how it says so
    instead of reporting the two as separate findings.
    """
    lost = {g: {m for m, v in statuses.get(g, {}).items() if v} for g in genes}
    out: list[dict[str, Any]] = []
    for i, a in enumerate(genes):
        for b in genes[i + 1 :]:
            sa, sb = lost.get(a) or set(), lost.get(b) or set()
            if not sa or not sb:
                continue
            overlap = len(sa & sb) / min(len(sa), len(sb))
            if overlap >= threshold:
                out.append(
                    {
                        "genes": [a, b],
                        "lost_lines": [len(sa), len(sb)],
                        "shared_lines": len(sa & sb),
                        "overlap_of_the_smaller_set": round(overlap, 3),
                        "reading": f"{a} and {b} are lost in the same lines; a pair test cannot say "
                        "which of them the dependency belongs to",
                    }
                )
    return sorted(out, key=lambda r: -r["overlap_of_the_smaller_set"])


def shuffled_null(
    pairs: list[tuple[str, str]],
    statuses: dict[str, dict[str, bool]],
    store: _RankStore,
    permutations: int = NULL_PERMUTATIONS,
    fdr: float = FDR,
    seed: int = SEED,
) -> dict[str, Any]:
    """The sweep again with the loss labels permuted: how many hits appear by chance.

    Group sizes are preserved per loss gene, so the null keeps the sweep's shape
    and destroys only the association. If a permuted sweep produces hits at the
    same rate as the real one, the real hits are the multiple-testing artefact
    they look like.
    """
    rng = random.Random(seed)
    coverages = {gene: store.canonical(frozenset(status)) for gene, status in statuses.items()}
    counts: list[int] = []
    for _ in range(permutations):
        permuted: dict[str, dict[str, bool]] = {}
        for gene, status in statuses.items():
            models = list(status)
            n_lost = sum(1 for v in status.values() if v)
            lost = set(rng.sample(models, n_lost)) if n_lost else set()
            permuted[gene] = {m: (m in lost) for m in models}
        ps: list[float] = []
        deltas: list[float] = []
        for lost_gene, dependency in pairs:
            status = permuted.get(lost_gene)
            if status is None:
                continue
            row = test_pair(lost_gene, dependency, "permuted", status, store, coverages[lost_gene])
            if "p" in row:
                ps.append(row["p"])
                deltas.append(row["cliffs_delta"])
        qs = benjamini_hochberg(ps)
        counts.append(sum(1 for q, d in zip(qs, deltas, strict=True) if q <= fdr and d < 0))
    return {
        "permutations": permutations,
        "seed": seed,
        "hits_per_permutation": counts,
        "mean_hits": round(sum(counts) / len(counts), 2) if counts else None,
        "max_hits": max(counts) if counts else None,
        "what_it_controls": "the sweep as a whole, not any one pair: group sizes are preserved and "
        "only the association is destroyed",
    }


def evidence_for(row: dict[str, Any]) -> Evidence:
    """A hit is a prediction, never a measurement of a therapy.

    Confidence is bounded at 0.7 whatever the statistics say, because the
    observation is a knockout in culture: strength of evidence for a target,
    not for an effect in a patient.
    """
    q = row.get("q")
    size = min(1.0, abs(row.get("cliffs_delta") or 0.0) / 0.5)
    significance = 0.0 if q is None else max(0.0, min(1.0, 1.0 - q / FDR))
    confidence = 0.2 + 0.5 * significance * size
    return prediction(
        f"GenomeOS synthetic-lethality sweep over {depmap.RELEASE}",
        f"{row['dependency']} dependency is stronger in cell lines with {row['lost']} "
        f"{row['loss'].replace('_', ' ')} (n={row.get('lines_lost')} lost vs {row.get('lines_intact')} "
        f"intact, Cliff's delta {row.get('cliffs_delta')}, q={None if q is None else round(q, 4)})",
        round(confidence, 3),
    )


def run(
    net: bool = True,
    log=None,
    cache_dir: Path = depmap.CACHE,
    permutations: int = NULL_PERMUTATIONS,
    fdr: float = FDR,
) -> dict[str, Any]:
    """The whole sweep: candidate space, controls, hits, and the permuted null."""
    paralogues = paralogue_table(list(PANEL), net=net, log=log)
    pairs, space = candidate_pairs(paralogues)
    dependencies = sorted({b for _a, b in pairs})
    loss_genes = sorted({a for a, _b in pairs} | {c["lost"] for c in PREREGISTERED if c["lost"] != "MSI"})
    if log:
        print(
            f"  {len(pairs)} candidate pairs: {len(loss_genes)} lost genes x {len(dependencies)} "
            f"dependencies",
            file=log,
            flush=True,
        )
    effect_table = depmap.matrix("gene_effect", set(dependencies), net=net, cache_dir=cache_dir, log=log)
    damaging_table = depmap.matrix("damaging", set(loss_genes), net=net, cache_dir=cache_dir, log=log)
    cn_table = depmap.matrix("absolute_cn", set(loss_genes), net=net, cache_dir=cache_dir, log=log)
    sig = depmap.signatures(net=net, cache_dir=cache_dir, log=log)
    model_rows = depmap.models(net=net, cache_dir=cache_dir, log=log)
    effect = effect_table["values"]
    damaging = damaging_table["values"]
    absolute_cn = cn_table["values"]
    msi = {m: v for m, v in ((m, _float(row.get("MSIScore"))) for m, row in sig.items()) if v is not None}
    store = _RankStore(effect)

    statuses: dict[str, dict[str, bool]] = {}
    for gene in loss_genes:
        status = loss_status(gene, "loss", damaging, absolute_cn, msi)
        if sum(1 for v in status.values() if v) >= MIN_LOST:
            statuses[gene] = status
    coverages = {gene: store.canonical(frozenset(status)) for gene, status in statuses.items()}

    rows: list[dict[str, Any]] = []
    for lost_gene, dependency in pairs:
        status = statuses.get(lost_gene)
        if status is None:
            continue
        rows.append(test_pair(lost_gene, dependency, "loss", status, store, coverages[lost_gene]))
    tested = [r for r in rows if "p" in r]
    qs = benjamini_hochberg([r["p"] for r in tested])
    for r, q in zip(tested, qs, strict=True):
        r["q"] = q
        r["verdict"] = (
            "recovered"
            if q <= fdr and r["cliffs_delta"] < 0
            else ("nominal" if r["p"] <= 0.05 and r["cliffs_delta"] < 0 else "not_recovered")
        )
    hits = sorted(
        (r for r in tested if r["verdict"] == "recovered"),
        key=lambda r: (r["q"], r["cliffs_delta"]),
    )

    controls = control_rows(damaging, absolute_cn, msi, store, tested, fdr)
    same_event = co_loss(statuses, sorted({r["lost"] for r in hits}))
    null = shuffled_null(
        [(a, b) for a, b in pairs if a in statuses],
        statuses,
        store,
        permutations=permutations,
        fdr=fdr,
    )
    lineages: dict[str, int] = {}
    for m in effect:
        lineage = (model_rows.get(m) or {}).get("OncotreeLineage") or "unknown"
        lineages[lineage] = lineages.get(lineage, 0) + 1
    return {
        "claim": CLAIM,
        "preregistered": [dict(c) for c in PREREGISTERED],
        "thresholds": {
            "fdr": fdr,
            "min_lines_lost": MIN_LOST,
            "min_lines_intact": MIN_INTACT,
            "damaging_min": DAMAGING_MIN,
            "deletion_max_copies": DELETION_MAX_COPIES,
            "msi_high_msisensor2": MSI_HIGH,
            "test": "Mann-Whitney rank sum, one-sided, normal approximation with tie and "
            "continuity correction",
            "multiple_testing": "Benjamini-Hochberg across every tested pair in the sweep",
            "verdict_rules": VERDICTS,
        },
        "candidate_space": space,
        "data": {
            "provenance": depmap.provenance(),
            "cell_lines_with_gene_effect": len(effect),
            "cell_lines_with_mutation_calls": len(damaging),
            "cell_lines_with_absolute_cn": len(absolute_cn),
            "cell_lines_with_msi_score": len(msi),
            "msi_high_lines": sum(1 for v in msi.values() if v >= MSI_HIGH),
            "dependency_genes_found": len(effect_table["genes"]),
            "dependency_genes_missing": effect_table["missing"][:20],
            "lineages": dict(sorted(lineages.items(), key=lambda kv: -kv[1])),
        },
        "tests": len(tested),
        "untestable": sum(1 for r in rows if r.get("verdict") == "untestable"),
        "loss_genes_with_enough_lines": len(statuses),
        "hits": len(hits),
        "controls": controls,
        "controls_recovered": sum(1 for c in controls if c["agrees_with_preregistration"]),
        "shuffled_null": null,
        "loss_genes_that_are_the_same_event": same_event,
        "confounders": [
            "a co-deleted neighbour is indistinguishable from the gene that carries the effect; "
            "loss_genes_that_are_the_same_event lists the pairs this sweep cannot separate",
            "a loss enriched in one lineage brings that lineage's dependencies with it; no lineage "
            "covariate is fitted here, so a hit can be a lineage effect wearing a gene's name",
            "significance is not effect size: with about a thousand lines a shift of 0.05 gene-effect "
            "units reaches p = 0.015, so every row carries the two group medians beside the p-value",
            "a damaging-mutation call is a prediction about a variant's consequence, not a measured "
            "absence of protein, and a heterozygous damaging mutation is counted as loss",
        ],
        "top_hits": [
            {**r, "q": round(r["q"], 5), "p": round(r["p"], 6), "evidence": evidence_for(r).to_dict()}
            for r in hits[:40]
        ],
        "not_a_therapy": (
            "A dependency is a hypothesis about a target. No molecule, construct, formulation or "
            "protocol follows from it here, and a cell-line dependency is not a clinical effect."
        ),
    }


def control_rows(
    damaging: dict[str, dict[str, float]],
    absolute_cn: dict[str, dict[str, float]],
    msi: dict[str, float],
    store: _RankStore,
    tested: list[dict[str, Any]],
    fdr: float = FDR,
) -> list[dict[str, Any]]:
    """The pre-registered controls, each under its own loss definition.

    The controls are a family declared in advance, so their p-values get their
    own Benjamini-Hochberg correction across the seven of them rather than
    borrowing the sweep's. The sweep's own q for the same pair, under the
    sweep's union loss definition, is reported beside each one.
    """
    rows: list[dict[str, Any]] = []
    for control in PREREGISTERED:
        kind = control["loss"]
        gene = control["lost"]
        status = loss_status(gene, kind, damaging, absolute_cn, msi)
        own = test_pair(gene, control["dependency"], kind, status, store)
        in_sweep = next(
            (
                r
                for r in tested
                if r["lost"] == gene and r["dependency"] == control["dependency"] and r["loss"] == "loss"
            ),
            None,
        )
        rows.append(
            {
                **{k: control[k] for k in ("lost", "loss", "dependency", "expect", "why")},
                "result": own,
                "in_sweep": None
                if in_sweep is None
                else {
                    k: in_sweep.get(k)
                    for k in (
                        "lines_lost",
                        "lines_intact",
                        "cliffs_delta",
                        "median_difference",
                        "p",
                        "q",
                        "verdict",
                    )
                },
            }
        )
    testable = [r for r in rows if "p" in r["result"]]
    qs = benjamini_hochberg([r["result"]["p"] for r in testable])
    for r, q in zip(testable, qs, strict=True):
        r["result"]["q_within_controls"] = round(q, 6)
        r["result"]["p"] = round(r["result"]["p"], 8)
    for r in rows:
        res = r["result"]
        q = res.get("q_within_controls")
        delta = res.get("cliffs_delta")
        if q is None:
            r["verdict"] = "untestable"
        elif q <= fdr and (delta or 0) < 0:
            r["verdict"] = "recovered"
        elif res["p"] <= 0.05 and (delta or 0) < 0:
            r["verdict"] = "nominal"
        else:
            r["verdict"] = "not_recovered"
        r["agrees_with_preregistration"] = (
            r["verdict"] == "recovered" if r["expect"] == "recovered" else r["verdict"] != "recovered"
        )
    return rows


def _float(raw: Any) -> float | None:
    try:
        return float(raw)
    except (TypeError, ValueError):
        return None
