"""cBioPortal client (public REST API, no key, standard library only).

https://www.cbioportal.org/api — studies, cancer types, clinical attributes
and mutations. Used to distil cancer knowledge (which genes are mutated in
which cancers, how often, at which residues) into small committed results,
following the project's stream-distil-discard rule. The MCP server at
mcp.cbioportal.org needs authentication and an MCP client; the REST API is
open, reproducible and dependency-free, so GenomeOS uses REST.
"""

from __future__ import annotations

import json
import time
import urllib.request
from collections.abc import Iterable
from typing import Any

BASE = "https://www.cbioportal.org/api"

# A compact panel of well-established cancer genes (drivers, tumour suppressors,
# oncogenes) used for distillation and for ranking a sample's variants. Sources:
# the Cancer Gene Census / OncoKB consensus lists; representative, not exhaustive.
DRIVER_PANEL: tuple[str, ...] = (
    "TP53",
    "KRAS",
    "NRAS",
    "HRAS",
    "BRAF",
    "EGFR",
    "ERBB2",
    "PIK3CA",
    "PTEN",
    "AKT1",
    "MTOR",
    "TSC1",
    "TSC2",
    "APC",
    "CTNNB1",
    "SMAD4",
    "TGFBR2",
    "RB1",
    "CDKN2A",
    "CDK4",
    "CCND1",
    "MYC",
    "MYCN",
    "NF1",
    "NF2",
    "VHL",
    "IDH1",
    "IDH2",
    "TET2",
    "DNMT3A",
    "ASXL1",
    "EZH2",
    "KMT2D",
    "KMT2C",
    "ARID1A",
    "ARID2",
    "SMARCA4",
    "CREBBP",
    "EP300",
    "BRCA1",
    "BRCA2",
    "ATM",
    "ATR",
    "CHEK2",
    "PALB2",
    "MLH1",
    "MSH2",
    "MSH6",
    "PMS2",
    "POLE",
    "NOTCH1",
    "NOTCH2",
    "FBXW7",
    "KIT",
    "PDGFRA",
    "FLT3",
    "JAK2",
    "STAT3",
    "ALK",
    "RET",
    "ROS1",
    "MET",
    "FGFR1",
    "FGFR2",
    "FGFR3",
    "CDH1",
    "GATA3",
    "ESR1",
    "AR",
    "FOXA1",
    "SPOP",
    "KEAP1",
    "NFE2L2",
    "STK11",
    "RBM10",
    "KDM6A",
    "ERBB3",
    "MAP2K1",
    "GNAS",
    "GNAQ",
    "GNA11",
    "SF3B1",
    "U2AF1",
    "SRSF2",
    "NPM1",
    "CEBPA",
    "RUNX1",
    "WT1",
    "MEN1",
    "DAXX",
    "ATRX",
    "TERT",
    "TP63",
    "SOX2",
    "NKX2-1",
    "AXIN1",
    "AXIN2",
    "RNF43",
    "ZNRF3",
    "BAP1",
    "PBRM1",
    "SETD2",
    "KDM5C",
    "CIC",
    "FUBP1",
    "H3-3A",
    "PTCH1",
    "SMO",
    "SUFU",
    "GLI1",
)


class CBioPortal:
    def __init__(self, base: str = BASE, timeout: int = 120, sleep: float = 0.2) -> None:
        self.base = base
        self.timeout = timeout
        self.sleep = sleep

    # ---- transport ---------------------------------------------------------

    def _get(self, path: str, **params: Any) -> Any:
        q = "&".join(f"{k}={v}" for k, v in params.items() if v is not None)
        url = f"{self.base}{path}" + (f"?{q}" if q else "")
        req = urllib.request.Request(
            url, headers={"Accept": "application/json", "User-Agent": "GenomeOS/0.1"}
        )
        with urllib.request.urlopen(req, timeout=self.timeout) as r:  # noqa: S310
            return json.load(r)

    def _post(self, path: str, body: Any, **params: Any) -> Any:
        q = "&".join(f"{k}={v}" for k, v in params.items() if v is not None)
        url = f"{self.base}{path}" + (f"?{q}" if q else "")
        data = json.dumps(body).encode()
        req = urllib.request.Request(
            url,
            data=data,
            method="POST",
            headers={
                "Accept": "application/json",
                "Content-Type": "application/json",
                "User-Agent": "GenomeOS/0.1",
            },
        )
        with urllib.request.urlopen(req, timeout=self.timeout) as r:  # noqa: S310
            return json.load(r)

    # ---- catalogue -----------------------------------------------------------

    def studies(self) -> list[dict]:
        return self._get("/studies", pageSize=5000, projection="SUMMARY")

    def cancer_types(self) -> list[dict]:
        return self._get("/cancer-types", pageSize=5000)

    def genes(self, symbols: Iterable[str]) -> dict[str, int]:
        """HGNC symbol -> Entrez id."""
        out = self._post("/genes/fetch", list(symbols), geneIdType="HUGO_GENE_SYMBOL", projection="SUMMARY")
        return {g["hugoGeneSymbol"]: g["entrezGeneId"] for g in out}

    def sample_attribute(self, study: str, attribute: str) -> dict[str, str]:
        """sampleId -> value of a clinical sample attribute (e.g. CANCER_TYPE)."""
        rows = self._get(
            f"/studies/{study}/clinical-data",
            clinicalDataType="SAMPLE",
            attributeId=attribute,
            pageSize=200000,
        )
        return {r["sampleId"]: r["value"] for r in rows}

    def mutations(self, profile: str, sample_list: str, entrez_ids: list[int]) -> list[dict]:
        """All mutation rows for the genes in a molecular profile / sample list."""
        return self._post(
            f"/molecular-profiles/{profile}/mutations/fetch",
            {"sampleListId": sample_list, "entrezGeneIds": entrez_ids},
            projection="SUMMARY",
            pageSize=10_000_000,
        )

    # ---- expression -----------------------------------------------------------

    #: Expression profiles, most informative first. The reference-normal z-score
    #: is the one worth having: it is tumour measured against matched normal
    #: tissue, so it is unit-free and says what a raw RSEM value cannot.
    EXPRESSION_PROFILES: tuple[tuple[str, str, str], ...] = (
        (
            "rna_seq_v2_mrna_median_all_sample_ref_normal_Zscores",
            "z_vs_normal",
            "z-score of each tumour against the study's normal samples",
        ),
        (
            "rna_seq_v2_mrna_median_Zscores",
            "z_vs_tumours",
            "z-score of each tumour against the other tumours in the study",
        ),
        ("rna_seq_v2_mrna", "rsem", "RSEM expression value, not comparable across datasets"),
    )

    def molecular_profiles(self, study: str) -> list[dict]:
        return self._get(f"/studies/{study}/molecular-profiles")

    def expression_profile(self, study: str) -> tuple[str, str, str] | None:
        """The best available expression profile of a study, and what it means."""
        available = {p["molecularProfileId"] for p in self.molecular_profiles(study)}
        for suffix, kind, meaning in self.EXPRESSION_PROFILES:
            if f"{study}_{suffix}" in available:
                return f"{study}_{suffix}", kind, meaning
        return None

    def expression(self, profile: str, sample_list: str, entrez_ids: list[int]) -> list[dict]:
        return self._post(
            f"/molecular-profiles/{profile}/molecular-data/fetch",
            {"sampleListId": sample_list, "entrezGeneIds": entrez_ids},
            projection="SUMMARY",
        )

    def expression_distribution(
        self, study: str, symbols: Iterable[str], batch: int = 60, progress=None
    ) -> dict[str, Any]:
        """Per gene, the distribution of its expression across a cohort's tumours.

        A cohort cannot say what one patient's tumour does. It says how often a
        gene is raised above normal tissue in this cancer type at all, which is
        the difference between a target worth investigating and a gene that
        happens to be mutated.
        """
        chosen = self.expression_profile(study)
        if chosen is None:
            raise ValueError(f"{study} has no mRNA expression profile")
        profile, kind, meaning = chosen
        sample_list = f"{study}_all"
        ids = self.genes(symbols)
        rev = {v: k for k, v in ids.items()}
        entrez = list(ids.values())
        values: dict[str, list[float]] = {}
        for i in range(0, len(entrez), batch):
            chunk = entrez[i : i + batch]
            for r in self.expression(profile, sample_list, chunk):
                v = r.get("value")
                sym = rev.get(r["entrezGeneId"])
                if v is None or sym is None:
                    continue
                values.setdefault(sym, []).append(float(v))
            if progress:
                progress(min(i + batch, len(entrez)), len(entrez))
            time.sleep(self.sleep)
        return {
            "study": study,
            "profile": profile,
            "kind": kind,
            "meaning": meaning,
            "samples": max((len(v) for v in values.values()), default=0),
            "genes": {g: summarise(v, kind) for g, v in sorted(values.items())},
            "missing": sorted(set(ids) - set(values)),
            "evidence": {
                "kind": "experimental",
                "source": f"cBioPortal {study} ({profile}), via public REST API",
            },
            "confidence": 0.8,
        }

    # ---- distillation ----------------------------------------------------------

    def distil_study(
        self, study: str, panel: Iterable[str] = DRIVER_PANEL, batch: int = 10, progress=None
    ) -> dict[str, Any]:
        """Per cancer type and per gene: fraction of samples mutated, top protein changes."""
        profile = f"{study}_mutations"
        sample_list = f"{study}_all"
        types = self.sample_attribute(study, "CANCER_TYPE")
        n_by_type: dict[str, int] = {}
        for t in types.values():
            n_by_type[t] = n_by_type.get(t, 0) + 1
        ids = self.genes(panel)
        symbols = {v: k for k, v in ids.items()}
        entrez = list(ids.values())
        per_gene: dict[str, dict[str, Any]] = {}
        for i in range(0, len(entrez), batch):
            chunk = entrez[i : i + batch]
            rows = self.mutations(profile, sample_list, chunk)
            for m in rows:
                sym = symbols.get(m["entrezGeneId"])
                if not sym:
                    continue
                g = per_gene.setdefault(sym, {"samples": set(), "by_type": {}, "changes": {}})
                sid = m["sampleId"]
                g["samples"].add(sid)
                ct = types.get(sid, "Unknown")
                g["by_type"].setdefault(ct, set()).add(sid)
                pc = m.get("proteinChange") or "?"
                g["changes"][pc] = g["changes"].get(pc, 0) + 1
            if progress:
                progress(min(i + batch, len(entrez)), len(entrez))
            time.sleep(self.sleep)
        n = len(types)
        out_genes = {}
        for sym, g in per_gene.items():
            out_genes[sym] = {
                "samples_mutated": len(g["samples"]),
                "frequency": round(len(g["samples"]) / n, 4),
                "by_cancer_type": {
                    ct: round(len(s) / n_by_type[ct], 4)
                    for ct, s in g["by_type"].items()
                    if n_by_type.get(ct, 0) >= 20
                },
                "hotspots": sorted(g["changes"].items(), key=lambda kv: -kv[1])[:8],
            }
        return {
            "study": study,
            "samples": n,
            "cancer_types": {k: v for k, v in sorted(n_by_type.items(), key=lambda kv: -kv[1])},
            "panel": list(ids),
            "genes": out_genes,
            "evidence": {
                "kind": "experimental",
                "source": f"cBioPortal {study} (tumour sequencing), via public REST API",
            },
            "confidence": 0.8,
        }


def summarise(values: list[float], kind: str, raised: float = 2.0) -> dict[str, Any]:
    """The shape of one gene's expression across a cohort, kept to a few numbers.

    `fraction_raised` is only meaningful for a z-score profile, where 2.0 is the
    usual threshold for calling a tumour over-expressed; for a raw RSEM profile
    it is left out rather than computed from an arbitrary cutoff.
    """
    v = sorted(values)
    n = len(v)

    def q(f: float) -> float:
        if n == 1:
            return v[0]
        i = f * (n - 1)
        lo = int(i)
        hi = min(lo + 1, n - 1)
        return v[lo] + (v[hi] - v[lo]) * (i - lo)

    out: dict[str, Any] = {
        "n": n,
        "min": round(v[0], 3),
        "q1": round(q(0.25), 3),
        "median": round(q(0.5), 3),
        "q3": round(q(0.75), 3),
        "p90": round(q(0.90), 3),
        "max": round(v[-1], 3),
    }
    if kind.startswith("z_"):
        out["fraction_raised"] = round(sum(1 for x in v if x >= raised) / n, 4)
        out["raised_threshold"] = raised
    return out


def percentile_of(value: float, summary: dict[str, Any]) -> float | None:
    """Roughly where a value falls in a distribution, from its quantiles only.

    Interpolates between the stored quantiles; it is a position, not a test,
    and it is only valid when the value is on the same scale as the cohort.
    """
    points = [
        (summary.get("min"), 0.0),
        (summary.get("q1"), 0.25),
        (summary.get("median"), 0.5),
        (summary.get("q3"), 0.75),
        (summary.get("p90"), 0.90),
        (summary.get("max"), 1.0),
    ]
    points = [(x, f) for x, f in points if x is not None]
    if len(points) < 2:
        return None
    if value <= points[0][0]:
        return 0.0
    for (x0, f0), (x1, f1) in zip(points, points[1:], strict=False):
        if value <= x1:
            if x1 == x0:
                return round(f1, 3)
            return round(f0 + (f1 - f0) * (value - x0) / (x1 - x0), 3)
    return 1.0
