"""Epigenetic clocks (task 3.1) from published coefficients, no dependencies.

Horvath 2013 (353 CpGs, multi-tissue) and Hannum 2013 (71 CpGs, blood) are
elastic-net linear models on methylation beta values. Coefficient tables are
the ones redistributed by the biolearn project (data/knowledge/Horvath1.csv,
Hannum.csv). Horvath's model predicts a transformed age:

    F(age) = log(age + 1) - log(adult + 1)         for age <= adult
           = (age - adult) / (adult + 1)           for age >  adult,  adult = 20

so the prediction is inverted with F^-1. Missing CpGs are imputed with 0.5 and
reported, because a clock evaluated on half its sites is a weaker measurement.
"""

from __future__ import annotations

import csv
import gzip
import math
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path

from genomeos.ir import Evidence, EvidenceKind, Parameter

ADULT = 20.0
HORVATH_INTERCEPT = 0.695507258  # "(Intercept)" row of Horvath 2013 Additional file 3
CLOCK_EVIDENCE = {
    "horvath2013": Evidence(
        EvidenceKind.EXPERIMENTAL, "Horvath 2013, Genome Biology 14:R115", note="353 CpGs, multi-tissue"
    ),
    "hannum2013": Evidence(
        EvidenceKind.EXPERIMENTAL, "Hannum et al. 2013, Molecular Cell 49:359", note="71 CpGs, whole blood"
    ),
}


@dataclass(slots=True)
class Clock:
    name: str
    intercept: float
    coefficients: dict[str, float]
    transform: str = "linear"  # "horvath" applies the age transform

    @classmethod
    def from_csv(cls, path: str | Path, name: str, transform: str = "linear") -> Clock:
        intercept = 0.0
        coefs: dict[str, float] = {}
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh):
                cpg = (row.get("CpGmarker") or row.get("cpg") or row.get("Marker") or "").strip()
                val = float(
                    row.get("CoefficientTraining") or row.get("Coefficient") or row.get("coefficient") or 0.0
                )
                if cpg.lower() in ("(intercept)", "intercept"):
                    intercept = val
                elif cpg:
                    coefs[cpg] = val
        return cls(name, intercept, coefs, transform)

    @staticmethod
    def horvath(path: str | Path = "data/knowledge/Horvath1.csv") -> Clock:
        c = Clock.from_csv(path, "horvath2013", "horvath")
        if c.intercept == 0.0:
            c.intercept = HORVATH_INTERCEPT
        return c

    @staticmethod
    def hannum(path: str | Path = "data/knowledge/Hannum.csv") -> Clock:
        return Clock.from_csv(path, "hannum2013", "linear")  # Hannum's model has no intercept

    def predict(self, betas: dict[str, float]) -> ClockResult:
        total = self.intercept
        missing = 0
        for cpg, coef in self.coefficients.items():
            b = betas.get(cpg)
            if b is None or (isinstance(b, float) and math.isnan(b)):
                b = 0.5
                missing += 1
            total += coef * b
        age = _inverse_horvath(total) if self.transform == "horvath" else total
        coverage = 1 - missing / len(self.coefficients)
        ev = CLOCK_EVIDENCE.get(self.name, Evidence(EvidenceKind.CURATED, self.name))
        conf = 0.8 * coverage
        return ClockResult(
            self.name,
            age,
            coverage,
            missing,
            Parameter(f"epigenetic_age_{self.name}", round(age, 2), "years", ev, conf),
        )


def _inverse_horvath(x: float) -> float:
    return (ADULT + 1) * math.exp(x) - 1 if x < 0 else (ADULT + 1) * x + ADULT


def horvath_transform(age: float) -> float:
    return math.log(age + 1) - math.log(ADULT + 1) if age <= ADULT else (age - ADULT) / (ADULT + 1)


@dataclass(slots=True)
class ClockResult:
    clock: str
    age: float
    coverage: float
    missing: int
    parameter: Parameter


# ---- methylation input -------------------------------------------------


@dataclass(slots=True)
class MethylationMatrix:
    samples: list[str]
    betas: dict[str, list[float]] = field(default_factory=dict)  # cpg -> per-sample values
    metadata: dict[str, dict[str, str]] = field(default_factory=dict)  # sample -> key -> value

    def sample(self, name: str) -> dict[str, float]:
        i = self.samples.index(name)
        return {cpg: vals[i] for cpg, vals in self.betas.items()}

    @classmethod
    def from_csv(cls, path: str | Path, delimiter: str = ",") -> MethylationMatrix:
        """First column CpG ids, remaining columns samples."""
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt") as fh:
            reader = csv.reader(fh, delimiter=delimiter)
            header = next(reader)
            m = cls(samples=[h.strip().strip('"') for h in header[1:]])
            for row in reader:
                if not row:
                    continue
                m.betas[row[0].strip().strip('"')] = [
                    float(x) if x not in ("", "NA", "null") else math.nan for x in row[1:]
                ]
        return m

    @classmethod
    def from_geo_series_matrix(cls, path: str | Path, cpgs: Iterable[str] | None = None) -> MethylationMatrix:
        """GEO series_matrix.txt(.gz): metadata lines start with '!', then a tab table."""
        wanted = set(cpgs) if cpgs is not None else None
        opener = gzip.open if str(path).endswith(".gz") else open
        meta: dict[str, list[str]] = {}
        m: MethylationMatrix | None = None
        with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("!"):
                    key, _, rest = line.rstrip("\n").partition("\t")
                    vals = [v.strip('"') for v in rest.split("\t")]
                    if key == "!Sample_geo_accession":
                        m = cls(samples=vals)
                    elif key == "!Sample_characteristics_ch1":
                        meta.setdefault("characteristics", []).append("|".join(vals))
                        for s, v in zip(m.samples if m else [], vals, strict=False):
                            k, _, val = v.partition(":")
                            m.metadata.setdefault(s, {})[k.strip().lower()] = val.strip()
                    continue
                if line.startswith('"ID_REF"') or line.startswith("ID_REF"):
                    continue
                if line.startswith("!series_matrix_table_end"):
                    break
                cpg, _, rest = line.rstrip("\n").partition("\t")
                cpg = cpg.strip('"')
                if wanted is not None and cpg not in wanted:
                    continue
                if m is None:
                    raise ValueError("no sample accession line before the table")
                m.betas[cpg] = [
                    float(x) if x not in ("", "null", "NA") else math.nan for x in rest.split("\t")
                ]
        if m is None:
            raise ValueError("not a GEO series matrix")
        return m


def pearson(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys, strict=True))
    return sxy / math.sqrt(sxx * syy) if sxx and syy else 0.0
