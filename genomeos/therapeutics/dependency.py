# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""What a deep deletion is worth: the dependency it creates, as evidence beside the targets.

A homozygously deleted gene is never a target: the tumour makes none of its
product (`cancer/alterations.py`, `GeneAlteration.removes_product`). What the
deletion points at is the gene the loss makes indispensable, and
`synthetic_lethality.py` measured those pairs across DepMap cell lines on
2026-09-16. This module is the wiring registered on 2026-10-10
(docs/THERAPEUTICS.md, "The dependency a deep deletion creates"):

* the trigger is a `deep_deletion` call in this patient's copy-number table,
  and nothing else;
* the partners are the sweep's 44 recovered pairs, frozen below, and no others;
* each partner carries the registered sweep row and a new lineage test, the
  DepMap lines of the lineage the caller states, deleted against kept;
* the result is a separate list, `dependency_candidates`, which the target
  ranking never reads, so no target candidate can move because of it.

A cell-line dependency is a hypothesis for a tumour, never a measurement in it
and never evidence that a therapy works. No molecule, construct, dose or
protocol follows from anything here.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any

from genomeos.cancer.alterations import DEEP_LOSS_COPIES, GeneAlteration

from . import depmap
from . import synthetic_lethality as sl
from .evidence import Evidence, today
from .evidence import patient as patient_evidence

#: Where the partners come from: the sweep's commit and its committed result.
SWEEP_COMMIT = "2305bee"
SWEEP_RESULT = "data/results/synthetic_lethality.json"

#: The 44 pairs the sweep recovered at q <= 0.10, in the sweep's own order:
#: (lost gene, dependency, lines lost, lines kept, median gene effect in the lost
#: lines, median in the kept lines, Cliff's delta, BH q across 19,749 tests).
#: The first 40 are the committed result's `top_hits`; the last four were
#: reproduced by re-running the sweep offline from the DepMap cache, which gave
#: back the committed 19,749 tests and 44 hits. Loss is the sweep's union
#: definition (a damaging mutation or fewer than 0.5 copies), pan-cancer.
REGISTERED_PAIRS: tuple[tuple[str, str, int, int, float, float, float, float], ...] = (
    ("MTAP", "WDR77", 177, 771, -1.635, -1.288, -0.441, 5.554e-16),
    ("ARID1A", "ARID1B", 121, 865, -0.116, 0.021, -0.428, 1.191e-10),
    ("MTAP", "PRMT5", 177, 771, -1.132, -0.975, -0.36, 2.643e-10),
    ("CDKN2B", "WDR77", 274, 675, -1.472, -1.279, -0.286, 1.134e-08),
    ("CDKN2A", "WDR77", 427, 538, -1.417, -1.272, -0.245, 1.163e-07),
    ("ARID1A", "WRN", 121, 865, -0.136, -0.049, -0.339, 2.511e-06),
    ("SMARCA4", "SMARCA2", 63, 899, -0.122, 0.013, -0.443, 5.734e-06),
    ("EP300", "CREBBP", 81, 884, -0.188, 0.028, -0.393, 5.771e-06),
    ("KMT2D", "RAD50", 103, 880, -0.122, 0.006, -0.324, 7.91e-05),
    ("KMT2D", "WRN", 103, 880, -0.144, -0.05, -0.323, 7.91e-05),
    ("ARID1A", "ARID2", 121, 865, -0.458, -0.337, -0.286, 0.0003111),
    ("MSH3", "WRN", 32, 925, -0.311, -0.049, -0.496, 0.001461),
    ("ARID1A", "KDM6A", 121, 865, 0.045, 0.106, -0.26, 0.002597),
    ("RECQL", "WRN", 10, 941, -0.901, -0.05, -0.823, 0.005202),
    ("SMARCA4", "RBBP8", 63, 899, -1.028, -0.901, -0.331, 0.006709),
    ("CDKN2B", "PRMT5", 274, 675, -1.049, -0.992, -0.182, 0.006709),
    ("KMT2D", "RIF1", 103, 880, -0.396, -0.325, -0.264, 0.006733),
    ("CDKN2A", "PRMT5", 427, 538, -1.044, -0.974, -0.163, 0.007705),
    ("RAD50", "WRN", 27, 929, -0.783, -0.05, -0.482, 0.01007),
    ("CREBBP", "EP300", 56, 901, -0.489, -0.253, -0.335, 0.01228),
    ("TP53BP1", "RBBP8", 28, 929, -1.088, -0.906, -0.459, 0.01601),
    ("MSH6", "KDM6A", 42, 917, 0.016, 0.103, -0.367, 0.02558),
    ("ATAD5", "WRN", 16, 938, -0.498, -0.05, -0.576, 0.03303),
    ("ARID1B", "WRN", 42, 923, -0.174, -0.05, -0.354, 0.03832),
    ("KMT2C", "KMT2D", 58, 911, -0.473, -0.331, -0.305, 0.03832),
    ("ARID1A", "RBBP8", 121, 865, -0.983, -0.9, -0.218, 0.03832),
    ("MTAP", "MAT2A", 177, 771, -1.867, -1.724, -0.186, 0.03924),
    ("RAD51AP1", "WRN", 11, 940, -0.433, -0.05, -0.669, 0.0456),
    ("BRCA1", "HELQ", 21, 934, -0.148, -0.061, -0.487, 0.0456),
    ("EXO1", "PARP2", 9, 944, -0.2, -0.077, -0.731, 0.04901),
    ("MSH3", "XRCC4", 32, 925, -0.352, -0.261, -0.392, 0.04901),
    ("KMT2D", "MRE11", 103, 880, -0.461, -0.355, -0.228, 0.04901),
    ("BRCA1", "RAD54B", 21, 934, -0.137, -0.085, -0.479, 0.05152),
    ("POLQ", "CDKN2B", 19, 923, -0.089, 0.044, -0.498, 0.05777),
    ("USP1", "ZSWIM7", 7, 945, -0.251, -0.045, -0.812, 0.05835),
    ("BRCA2", "WRN", 37, 923, -0.175, -0.05, -0.359, 0.05835),
    ("MSH3", "PNKP", 32, 925, -0.517, -0.381, -0.381, 0.06633),
    ("FANCA", "RAD50", 16, 937, -0.175, 0.004, -0.53, 0.07132),
    ("ERCC5", "SLFN11", 6, 945, -0.21, -0.021, -0.852, 0.07987),
    ("BRCA1", "CENPX", 21, 934, -0.38, -0.224, -0.455, 0.08762),
    ("DNMT1", "RAD50", 10, 938, -0.186, 0.004, -0.652, 0.0903),
    ("SMARCA4", "MAD2L2", 63, 899, -1.026, -0.818, -0.267, 0.0903),
    ("BRCA1", "RAD23A", 21, 934, -0.282, -0.188, -0.45, 0.09474),
    ("ATM", "POLH", 30, 930, -0.313, -0.189, -0.377, 0.09822),
)

#: The loss genes the sweep could test: at least 5 lines lost each.
_TESTABLE = """
    ARID1A ARID1B ARID2 ATAD1 ATAD5 ATM ATR ATRIP BARD1 BLM BRCA1 BRCA2 BRD7 BRIP1 CABP2 CALM2
    CDC25B CDC25C CDC7 CDKN2A CDKN2B CETN3 CHD1 CHD1L CHD2 CHEK1 CHEK2 CLSPN CREBBP CSNK1A1L
    CSNK1D DBF4 DBF4B DCLRE1C DMC1 DNA2 DNMT1 DNMT3A DNMT3B DNTT DPF3 EED EME2 EP300 ERCC2
    ERCC3 ERCC4 ERCC5 ERCC6 EXO1 EZH1 EZH2 FANCA FANCC FANCD2 FANCG FANCI FANCL FANCM FIGNL1
    GEN1 GTF2H1 HELLS HELQ KAT6A KAT6B KATNAL2 KDM6A KDM6B KMT2C KMT2D LIG3 LIG4 MAD2L2 MAT1A
    MLH1 MLH3 MRE11 MSH2 MSH3 MSH4 MSH6 MTAP MTREX MUS81 NBN NEIL3 NHEJ1 PALB2 PARP1 PBRM1
    PIF1 PMS1 PMS2 PNKP POLA1 POLB POLD1 POLH POLI POLK POLN POLQ PRKDC RAD18 RAD23B RAD50
    RAD51AP1 RAD51B RAD54B RAD9B RBBP8 RBBP8NL RECQL RECQL4 RECQL5 REV1 REV3L RFWD3 RIF1 RIOK3
    RMI1 RMI2 RPA1 RPA4 RTEL1 SAMHD1 SETX SLFN11 SLFN12L SLFN13 SLFN14 SLFN5 SLX4 SMARCA1
    SMARCA2 SMARCA4 SMARCAD1 SMARCB1 SMARCC1 SMARCC2 SPAST SPRTN STAG1 STAG2 STAG3 SUZ12 TDG
    TDP1 TIMELESS TOP1 TOP2A TOP2B TOP3A TOP3B TOPBP1 TP53BP1 TRDMT1 TTF2 UBE2I USP1 WDR48
    WEE2 WRN WRNIP1 XPC XRCC2 XRCC5 XRCC6 ZNF451 ZRANB3
    """
SWEEP_TESTABLE: frozenset[str] = frozenset(_TESTABLE.split())

#: The loss genes the sweep asked about and could not test: fewer than 5 lines lost each.
_ASKED_UNTESTABLE = """
    APEX1 APEX2 APLF APTX ASF1A ASF1B BRD9 CABP5 CABP7 CALM1 CALM3 CALML3 CALML4 CALML5 CALML6
    CDC25A CDC34 CENPS CENPX CETN1 CETN2 CSNK1A1 DCLRE1B DDB1 DDB2 DDX11 DNMT3L DPF1 DPF2 DRC8
    EME1 ERCC1 ERCC8 FANCB FANCE FANCF FEN1 GTF2H5 HINT3 HLTF HMG20B HMGB3 HUS1 HUS1B KAT5
    KAT7 KAT8 KATNA1 KATNAL1 LIG1 MAD2L1 MAT2A MBD4 MECP2 MPG MSH5 MUTYH NEIL1 NEIL2 NTHL1
    NUCKS1 OGG1 PARG PARP2 PARP3 PCNA PDIK1L PKMYT1 PNP POLL POLM PRMT5 RAD1 RAD17 RAD23A
    RAD51 RAD51C RAD51D RAD52 RAD54L RAD9A RIOK1 RNASEH2A RNASEH2B RNASEH2C RPA2 RPA3 SHLD1
    SHLD2 SHLD3 SLFN12 SMARCA5 SMARCAL1 SMARCE1 SMUG1 SNAI3 SUV39H1 SUV39H2 SWSAP1 TDP2 TIPIN
    TOP1MT UBE2A UBE2B UBE2C UBE2F UBE2G1 UBE2G2 UBE2H UBE2K UBE2N UBE2NL UBE2T UBE2V1 UBE2V2
    UBE2W UNG USP18 USP27X USP50 UTY VPS4A VPS4B WDR77 WEE1 XNDC1N XPA XRCC1 XRCC3 XRCC4 YBEY
    ZNF575 ZSWIM7
    """
SWEEP_ASKED_UNTESTABLE: frozenset[str] = frozenset(_ASKED_UNTESTABLE.split())

SWEEP_ASKED: frozenset[str] = SWEEP_TESTABLE | SWEEP_ASKED_UNTESTABLE

#: Loss genes the sweep found lost in the same lines (`loss_genes_that_are_the_same_event`):
#: a pair test cannot say which of the two carries the dependency.
SAME_EVENT: tuple[tuple[str, str, float], ...] = (
    ("CDKN2A", "CDKN2B", 0.974),
    ("CDKN2A", "MTAP", 0.951),
    ("CDKN2B", "MTAP", 0.906),
)

STATUSES = {
    "partners_registered": "the sweep recovered at least one dependency for this loss at q <= 0.10",
    "asked_no_partner": "the sweep tested this loss and no pair reached q <= 0.10",
    "asked_untestable": "the sweep asked about this loss and fewer than 5 cell lines had lost the gene",
    "not_asked": "outside the sweep's candidate space: the sweep did not ask, which is not evidence "
    "that no dependency exists",
}

#: The lineage verdicts, in the order the list is sorted.
VERDICTS: dict[str, str] = {
    "supported_in_lineage": f"q <= {sl.FDR} across this analysis's lineage tests and the dependency is "
    "stronger in the deleted lines",
    "nominal_in_lineage": f"p <= 0.05 in the same direction, q above {sl.FDR}",
    "unknown": "the lineage test could not be run; the reason says why",
    "not_supported_in_lineage": "enough lines, and neither of the above",
    "contradicted_by_this_tumour": "this tumour deleted the partner too, and grows without it",
}
VERDICT_ORDER: tuple[str, ...] = tuple(VERDICTS)

CLOSING = (
    "A cell-line dependency is a hypothesis for this tumour, not a measurement in it, and not "
    "evidence that any therapy works."
)

ROUTE = (
    "a dependency partner would be reached by an inhibitor or a degrader, which GenomeOS does not "
    "model; it is listed with its evidence and never scored on the target dimensions"
)

LIMITS: tuple[str, ...] = (
    "the dependency is measured in cell lines grown in culture: no stroma, no immune system, and a "
    "growth-selected subset of the disease",
    "a CRISPR knockout removes the gene; an inhibitor or a degrader does not reproduce that, which is "
    "why BRCA1/2 against PARP1 reads weak in the same screens",
    "the sweep fitted no lineage covariate, so a pan-cancer pair can be a lineage effect; the lineage "
    "test is the check, where the lineage is stated and has the lines",
    "a deep-deletion call in this tumour comes from the supplied table and is not re-measured here",
)


@dataclass
class DepMapView:
    """The few DepMap columns the lineage test reads, from the local cache only."""

    effect: dict[str, dict[str, float]] = field(default_factory=dict)
    copies: dict[str, dict[str, float]] = field(default_factory=dict)
    copy_genes: frozenset[str] = frozenset()
    lineage_of: dict[str, str] = field(default_factory=dict)
    available: bool = False
    reason: str = ""

    @property
    def lineages(self) -> list[str]:
        return sorted({v for v in self.lineage_of.values() if v})


@lru_cache(maxsize=4)
def load_view(cache_dir: str = str(depmap.CACHE)) -> DepMapView:
    """Gene effect for the partners, copy number for the lost genes, lineage per line.

    Never streams: a pipeline run reads what the sweep cached and nothing else,
    so a missing cache is an `unknown` with its reason and not a download.
    """
    partners = {b for _a, b, *_ in REGISTERED_PAIRS}
    lost = {a for a, *_ in REGISTERED_PAIRS} | {g for a, b, _ in SAME_EVENT for g in (a, b)}
    try:
        effect = depmap.matrix("gene_effect", partners, net=False, cache_dir=Path(cache_dir))
        copies = depmap.matrix("absolute_cn", lost, net=False, cache_dir=Path(cache_dir))
        models = depmap.models(net=False, cache_dir=Path(cache_dir))
    except depmap.DepMapUnavailableError as e:
        return DepMapView(reason=f"the DepMap cache is not available on this machine ({e})")
    return DepMapView(
        effect=effect["values"],
        copies=copies["values"],
        copy_genes=frozenset(copies.get("genes") or ()),
        lineage_of={m: (row.get("OncotreeLineage") or "") for m, row in models.items()},
        available=True,
    )


def status_of(gene: str) -> str:
    if any(a == gene for a, *_ in REGISTERED_PAIRS):
        return "partners_registered"
    if gene in SWEEP_TESTABLE:
        return "asked_no_partner"
    if gene in SWEEP_ASKED_UNTESTABLE:
        return "asked_untestable"
    return "not_asked"


def same_event_genes(gene: str) -> list[tuple[str, float]]:
    out = [(b, o) for a, b, o in SAME_EVENT if a == gene] + [(a, o) for a, b, o in SAME_EVENT if b == gene]
    return sorted(out, key=lambda t: -t[1])


def copy_state(gene: str, alterations: Iterable[GeneAlteration], copy_number: dict[str, float]) -> str:
    """This tumour's copy state for one gene: deleted, kept, or unknown.

    Deleted is a deep-deletion call. Kept is any other call on the gene or a
    count above the deep-loss line. A gene the tables do not mention is unknown,
    never assumed diploid.
    """
    calls = [a for a in alterations if a.gene.upper() == gene]
    if any(a.kind == "deep_deletion" for a in calls):
        return "deleted"
    if calls:
        return "kept"
    cn = copy_number.get(gene)
    if cn is None:
        return "unknown"
    return "deleted" if cn <= DEEP_LOSS_COPIES else "kept"


def resolve_lineage(stated: str, view: DepMapView) -> tuple[str | None, str]:
    if not (stated or "").strip():
        return None, (
            "the tumour's lineage was not stated, so no lineage-matched cell lines were read; state it "
            "as a DepMap OncotreeLineage (for example --lineage 'CNS/Brain')"
        )
    if not view.available:
        return None, view.reason
    match = next((n for n in view.lineages if n.lower() == stated.strip().lower()), None)
    if match is None:
        return None, (
            f"'{stated}' is not a DepMap OncotreeLineage; the lineages are: {', '.join(view.lineages)}"
        )
    return match, ""


def lineage_test(
    lost: str,
    partner: str,
    lineage: str,
    view: DepMapView,
    restrict: tuple[tuple[str, str], ...] = (),
) -> dict[str, Any]:
    """Do the lines of this lineage that deeply deleted `lost` depend more on `partner`?

    `restrict` holds (same-event gene, this tumour's state of it) for every
    same-event gene this tumour's table states: the deleted group then keeps
    only lines that match this tumour on each of them, so a CDKN2A-deleted
    tumour that kept MTAP is compared with lines like it.
    """
    if lost not in view.copy_genes:
        return {
            "verdict": "unknown",
            "reason": f"DepMap carries no absolute copy number for {lost} (masked or absent)",
        }
    status: dict[str, bool] = {}
    for model, name in view.lineage_of.items():
        if name != lineage:
            continue
        cn = (view.copies.get(model) or {}).get(lost)
        if cn is None:
            continue
        deleted = cn < sl.DELETION_MAX_COPIES
        if deleted and not all(
            (other := (view.copies.get(model) or {}).get(gene)) is not None
            and (other < sl.DELETION_MAX_COPIES) == (state == "deleted")
            for gene, state in restrict
        ):
            continue
        status[model] = deleted
    row = sl.test_pair(lost, partner, "deletion", status, sl._RankStore(view.effect))
    if "p" not in row:
        return {
            **{k: row.get(k) for k in ("lines_lost", "lines_intact")},
            "verdict": "unknown",
            "reason": f"in the {lineage} lines: {row.get('reason')} (the test needs at least "
            f"{sl.MIN_LOST} deleted and {sl.MIN_INTACT} kept)",
        }
    return {k: row[k] for k in row if k not in ("lost", "loss", "dependency")}


def lineage_verdict(p: float, q: float, delta: float) -> str:
    """The sweep's verdict rule, read over this analysis's lineage tests."""
    if q <= sl.FDR and delta < 0:
        return "supported_in_lineage"
    if p <= 0.05 and delta < 0:
        return "nominal_in_lineage"
    return "not_supported_in_lineage"


def _sweep_row(pair: tuple[str, str, int, int, float, float, float, float]) -> dict[str, Any]:
    a, b, n_lost, n_kept, m_lost, m_kept, delta, q = pair
    return {
        "lost": a,
        "loss": "loss",
        "dependency": b,
        "lines_lost": n_lost,
        "lines_intact": n_kept,
        "median_effect_lost": m_lost,
        "median_effect_intact": m_kept,
        "cliffs_delta": delta,
        "q": q,
        "scope": "pan-cancer, the sweep's union loss definition (a damaging mutation or fewer than "
        f"{sl.DELETION_MAX_COPIES} copies)",
        "source": f"{SWEEP_RESULT} (commit {SWEEP_COMMIT})",
    }


def _lineage_evidence(a: str, b: str, lineage: str, t: dict[str, Any]) -> Evidence:
    """The lineage row as evidence: in vitro, never patient-derived, at most 0.7."""
    q = t.get("q")
    size = min(1.0, abs(t.get("cliffs_delta") or 0.0) / 0.5)
    significance = 0.0 if q is None else max(0.0, min(1.0, 1.0 - q / sl.FDR))
    confidence = 0.2 + 0.5 * significance * size
    stronger = t["verdict"] == "supported_in_lineage"
    claim = (
        f"In {depmap.RELEASE} {lineage} cell lines, {b} knockout "
        + ("slowed growth more" if stronger else "did not clearly slow growth more")
        + f" in the {t['lines_lost']} that deeply deleted {a} than in the {t['lines_intact']} that "
        f"kept it (median gene effect {t['median_effect_lost']} against {t['median_effect_intact']}, "
        f"Cliff's delta {t['cliffs_delta']}, p {t['p']:.2g}, q {q:.2g})"
    )
    return Evidence(
        source=f"{depmap.RELEASE} CRISPR gene effect and absolute copy number; GenomeOS lineage test",
        source_type="database",
        claim=claim,
        level="in_vitro",
        confidence=round(confidence, 3),
        identifier=depmap.DOI,
        retrieved_at=today(),
        version=depmap.RELEASE,
    )


def _reading(c: dict[str, Any]) -> str:
    a, b, t = c["lost_gene"], c["gene"], c["lineage"]
    verdict, name = t["verdict"], t.get("name") or "the stated"
    if verdict == "contradicted_by_this_tumour":
        head = f"This tumour deleted {b} too and grows: it is not dependent on {b}."
    elif verdict == "supported_in_lineage":
        head = (
            f"In the {t['lines_lost']} DepMap {name} lines that deeply deleted {a}, knocking out {b} "
            f"slowed growth more than in the {t['lines_intact']} that kept {a} (median gene effect "
            f"{t['median_effect_lost']} against {t['median_effect_intact']})."
        )
    elif verdict == "nominal_in_lineage":
        head = (
            f"In the {name} lines, {b} knockout shifted the same way where {a} was deleted at p "
            f"{t['p']:.2g}, not past the false-discovery rate across this analysis's lineage tests; "
            "unconfirmed in this lineage."
        )
    elif verdict == "not_supported_in_lineage":
        head = (
            f"In the {name} lines the dependency on {b} is not stronger where {a} was deleted; the "
            "pan-cancer pair stands as registered and is not borne out here."
        )
    else:
        head = f"Lineage-matched evidence unknown: {t['reason']}."
    s = c["sweep"]
    pan = (
        f" Pan-cancer, the sweep found {b} dependency stronger in the {s['lines_lost']} lines that "
        f"lost {a} (Cliff's delta {s['cliffs_delta']}, q {s['q']:.2g})."
    )
    return head + pan + " " + CLOSING


def empty() -> dict[str, Any]:
    return {"dependencies": [], "dependency_candidates": [], "dependency_lineage": None}


def assess(
    alterations: Iterable[GeneAlteration],
    copy_number: dict[str, float] | None = None,
    lineage: str = "",
    target_genes: Iterable[str] = (),
    source: str = "",
    view: DepMapView | None = None,
) -> dict[str, Any]:
    """The deep deletions of one tumour, and the dependencies each one registers.

    Returns `dependencies` (one record per deleted gene, with its status),
    `dependency_candidates` (one per registered pair whose lost gene this tumour
    deleted) and `dependency_lineage` (what lineage was read, or why none).
    Nothing here touches a target candidate.
    """
    alterations = list(alterations)
    copy_number = {k.upper(): v for k, v in (copy_number or {}).items()}
    deleted = sorted({a.gene.upper() for a in alterations if a.kind == "deep_deletion"})
    if not deleted:
        return empty()
    targets = {g.upper() for g in target_genes}
    calls = {a.gene.upper(): a for a in alterations if a.kind == "deep_deletion"}

    records: list[dict[str, Any]] = []
    for gene in deleted:
        status = status_of(gene)
        a = calls[gene]
        records.append(
            {
                "gene": gene,
                "event": a.label(),
                "status": status,
                "status_reason": STATUSES[status],
                "same_event": [
                    {
                        "gene": other,
                        "overlap_in_cell_lines": overlap,
                        "state_in_this_tumour": copy_state(other, alterations, copy_number),
                    }
                    for other, overlap in same_event_genes(gene)
                ],
                "partners": [b for x, b, *_ in REGISTERED_PAIRS if x == gene],
                "evidence": [
                    patient_evidence(
                        f"copy number ({source or 'supplied table'})",
                        f"{a.label()} called in this tumour: both copies of {gene} are gone",
                        0.85,
                    ).to_dict()
                ],
            }
        )

    wanted = [p for p in REGISTERED_PAIRS if p[0] in set(deleted)]
    if not wanted:
        return {"dependencies": records, "dependency_candidates": [], "dependency_lineage": None}
    if view is None:
        # The cache is read only when a lineage is stated: unstated, nothing in it is consulted.
        view = load_view() if (lineage or "").strip() else DepMapView()
    name, why_not = resolve_lineage(lineage, view)

    candidates: list[dict[str, Any]] = []
    for pair in wanted:
        a, b = pair[0], pair[1]
        record = next(r for r in records if r["gene"] == a)
        restrict = tuple(
            (s["gene"], s["state_in_this_tumour"])
            for s in record["same_event"]
            if s["state_in_this_tumour"] != "unknown"
        )
        if copy_state(b, alterations, copy_number) == "deleted":
            t = {
                "verdict": "contradicted_by_this_tumour",
                "reason": f"{b} is itself deep-deleted in this tumour",
            }
        elif name is None:
            t = {"verdict": "unknown", "reason": why_not}
        else:
            t = lineage_test(a, b, name, view, restrict)
        t["name"] = name
        t["restricted_to_lines_matching_this_tumour_on"] = (
            [{"gene": g, "state": st} for g, st in restrict] if name else []
        )
        unknown_same_event = [
            s["gene"] for s in record["same_event"] if s["state_in_this_tumour"] == "unknown"
        ]
        sweep = _sweep_row(pair)
        limitations = list(LIMITS)
        if record["same_event"]:
            limitations.insert(
                0,
                f"{a} is lost with {', '.join(s['gene'] for s in record['same_event'])} in the same cell "
                "lines, so the pair test cannot say which of them carries the dependency"
                + (
                    f"; this tumour's state of {', '.join(unknown_same_event)} is not in the supplied "
                    "table, so the cell-line comparison mixes lines that lost it and lines that kept it"
                    if unknown_same_event
                    else ""
                ),
            )
        candidates.append(
            {
                "gene": b,
                "lost_gene": a,
                "created_by": record["event"],
                "route": ROUTE,
                "also_a_target_candidate": b in targets,
                "partner_state_in_this_tumour": copy_state(b, alterations, copy_number),
                "sweep": sweep,
                "lineage": t,
                "evidence": [sl.evidence_for(sweep).to_dict()],
                "limitations": limitations,
            }
        )

    tested = [c for c in candidates if "p" in c["lineage"]]
    qs = sl.benjamini_hochberg([c["lineage"]["p"] for c in tested])
    for c, q in zip(tested, qs, strict=True):
        t = c["lineage"]
        t["q"] = q
        t["tests_in_this_analysis"] = len(tested)
        t["verdict"] = lineage_verdict(t["p"], q, t["cliffs_delta"])
        c["evidence"].append(_lineage_evidence(c["lost_gene"], c["gene"], t["name"], t).to_dict())
        t["p"] = round(t["p"], 8)
        t["q"] = round(q, 6)
    for c in candidates:
        c["lineage"]["verdict_rule"] = VERDICTS[c["lineage"]["verdict"]]
        c["reading"] = _reading(c)
    candidates.sort(
        key=lambda c: (
            VERDICT_ORDER.index(c["lineage"]["verdict"]),
            c["sweep"]["q"],
            c["sweep"]["cliffs_delta"],
            c["gene"],
        )
    )
    return {
        "dependencies": records,
        "dependency_candidates": candidates,
        "dependency_lineage": {
            "stated": lineage or None,
            "read": name,
            "reason": why_not or None,
            "release": depmap.RELEASE,
            "doi": depmap.DOI,
            "licence": depmap.LICENCE,
        },
    }
