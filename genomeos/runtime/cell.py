"""Cell ageing runtime: the "timer" in a cell.

A stochastic, per-cell state model with four ageing clocks:
    telomere length      shortens with each division (Hayflick limit)
    somatic mutations    accumulate roughly linearly with age
    epigenetic age       drifts with chronological age (methylation clocks)
    senescence           triggered by critically short telomeres or damage

Every numeric parameter is a BioIR Parameter with Evidence and Confidence.
The numbers below are order-of-magnitude values from the literature and are
marked accordingly; the point of v0.1 is the architecture (state separate
from genome, evidence attached to every rule), not precise prediction.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass

from genomeos.ir import Evidence, EvidenceKind, Parameter

E_HARLEY = Evidence(
    EvidenceKind.EXPERIMENTAL,
    "Harley, Futcher & Greider 1990, Nature 345:458",
    note="telomere loss per division in cultured fibroblasts",
)
E_HAYFLICK = Evidence(EvidenceKind.EXPERIMENTAL, "Hayflick & Moorhead 1961, Exp Cell Res 25:585")
E_CAGAN = Evidence(
    EvidenceKind.EXPERIMENTAL,
    "Cagan et al. 2022, Nature 604:517",
    note="somatic mutation rate per year across mammalian tissues",
)
E_HORVATH = Evidence(
    EvidenceKind.CURATED,
    "Horvath 2013, Genome Biology 14:R115",
    note="multi-tissue epigenetic clock; pace ~1 year/year in healthy adults",
)
E_TELOMERE_BIRTH = Evidence(
    EvidenceKind.EXPERIMENTAL,
    "Factor-Litvak et al. 2016, Pediatrics 137:e20153927",
    note="newborn leukocyte telomere length ~ 9-11 kb",
)
E_EST = Evidence(EvidenceKind.INFERRED, "GenomeOS v0.1 estimate", note="order of magnitude only")
E_HIYAMA = Evidence(
    EvidenceKind.CURATED,
    "Hiyama & Hiyama 2007, Br J Cancer 96:1020",
    note="telomerase activity in adult stem cells slows but does not stop attrition",
)
E_LTL = Evidence(
    EvidenceKind.EXPERIMENTAL,
    "Codd et al. 2021, Nature Genetics 53:1425",
    note="UK Biobank leukocyte telomere length declines ~20-30 bp/year",
)


@dataclass(slots=True)
class CellTypeParams:
    name: str
    divisions_per_year: Parameter
    telomere_loss_per_division_bp: Parameter
    mutations_per_year: Parameter
    telomerase_compensation: Parameter  # 0 = none, 1 = full maintenance


CELL_TYPES: dict[str, CellTypeParams] = {
    # Net attrition is calibrated so that blood/skin lose ~25-35 bp/year, the
    # range measured in large cohorts (E_LTL). Division rates are in-vivo
    # estimates (most adult fibroblasts are quiescent) and remain low confidence.
    "fibroblast": CellTypeParams(
        "fibroblast",
        Parameter("divisions_per_year", 0.5, "1/yr", E_EST, 0.3),
        Parameter("telomere_loss_per_division_bp", 70.0, "bp", E_HARLEY, 0.8),
        Parameter("mutations_per_year", 30.0, "1/yr", E_CAGAN, 0.6),
        Parameter("telomerase_compensation", 0.0, "fraction", E_HIYAMA, 0.7),
    ),
    "hematopoietic_stem": CellTypeParams(
        "hematopoietic_stem",
        Parameter("divisions_per_year", 1.3, "1/yr", E_EST, 0.3),
        Parameter("telomere_loss_per_division_bp", 60.0, "bp", E_HARLEY, 0.6),
        Parameter("mutations_per_year", 17.0, "1/yr", E_CAGAN, 0.6),
        Parameter("telomerase_compensation", 0.6, "fraction", E_HIYAMA, 0.4),
    ),
    "colon_crypt_stem": CellTypeParams(
        "colon_crypt_stem",
        Parameter("divisions_per_year", 60.0, "1/yr", E_EST, 0.4),
        Parameter("telomere_loss_per_division_bp", 30.0, "bp", E_EST, 0.3),
        Parameter("mutations_per_year", 43.0, "1/yr", E_CAGAN, 0.7),
        Parameter("telomerase_compensation", 0.98, "fraction", E_HIYAMA, 0.4),
    ),
    "neuron": CellTypeParams(
        "neuron",
        Parameter("divisions_per_year", 0.0, "1/yr", E_HAYFLICK, 0.9),
        Parameter("telomere_loss_per_division_bp", 0.0, "bp", E_EST, 0.5),
        Parameter("mutations_per_year", 20.0, "1/yr", E_CAGAN, 0.6),
        Parameter("telomerase_compensation", 0.0, "fraction", E_HIYAMA, 0.7),
    ),
}

GLOBAL_PARAMS: dict[str, Parameter] = {
    "telomere_at_birth_bp": Parameter("telomere_at_birth_bp", 10000.0, "bp", E_TELOMERE_BIRTH, 0.7),
    "telomere_sd_at_birth_bp": Parameter("telomere_sd_at_birth_bp", 700.0, "bp", E_EST, 0.3),
    # Senescence is a hazard, not a wall: the shortest of a cell's 92 telomere
    # ends triggers it long before the mean reaches the Hayflick range, so the
    # hazard rises exponentially as the mean approaches telomere_senescence_bp.
    "telomere_senescence_bp": Parameter("telomere_senescence_bp", 4000.0, "bp", E_HAYFLICK, 0.5),
    "telomere_hazard_scale_bp": Parameter("telomere_hazard_scale_bp", 1000.0, "bp", E_EST, 0.2),
    "senescence_hazard_base": Parameter("senescence_hazard_base", 0.02, "1/yr", E_EST, 0.2),
    "mutation_senescence_threshold": Parameter("mutation_senescence_threshold", 6000.0, "count", E_EST, 0.2),
    "epigenetic_pace_sd": Parameter("epigenetic_pace_sd", 0.08, "yr/yr", E_HORVATH, 0.5),
}


@dataclass(slots=True)
class Environment:
    """Exposures that modulate the clocks. 1.0 = population average."""

    proliferation_factor: float = 1.0  # e.g. chronic inflammation > 1
    mutagen_factor: float = 1.0  # e.g. smoking, UV > 1
    epigenetic_pace: float = 1.0  # biological years per chronological year


@dataclass(slots=True)
class CellState:
    cell_type: str
    age_years: float = 0.0
    divisions: int = 0
    telomere_bp: float = 10000.0
    somatic_mutations: int = 0
    epigenetic_age: float = 0.0
    senescent: bool = False
    pace: float = 1.0  # this cell's own epigenetic pace


def _poisson(rng: random.Random, lam: float) -> int:
    """Knuth's algorithm; fine for lam < ~500, which is all we need per tick."""
    if lam <= 0:
        return 0
    if lam > 500:  # normal approximation
        return max(0, int(round(rng.gauss(lam, math.sqrt(lam)))))
    limit = math.exp(-lam)
    k, p = 0, 1.0
    while True:
        k += 1
        p *= rng.random()
        if p <= limit:
            return k - 1


@dataclass(slots=True)
class TissueReport:
    age_years: float
    n_cells: int
    senescent_fraction: float
    mean_telomere_bp: float
    mean_mutations: float
    mean_epigenetic_age: float


class CellRuntime:
    def __init__(self, seed: int | None = None, env: Environment | None = None) -> None:
        self.rng = random.Random(seed)
        self.env = env or Environment()
        self.params = {k: p.value for k, p in GLOBAL_PARAMS.items()}

    def new_cell(self, cell_type: str) -> CellState:
        if cell_type not in CELL_TYPES:
            raise KeyError(f"unknown cell type {cell_type!r}; known: {sorted(CELL_TYPES)}")
        return CellState(
            cell_type=cell_type,
            telomere_bp=max(
                2000.0,
                self.rng.gauss(self.params["telomere_at_birth_bp"], self.params["telomere_sd_at_birth_bp"]),
            ),
            pace=max(0.3, self.rng.gauss(self.env.epigenetic_pace, self.params["epigenetic_pace_sd"])),
        )

    def tick(self, cell: CellState, dt_years: float) -> None:
        p = CELL_TYPES[cell.cell_type]
        cell.age_years += dt_years
        cell.epigenetic_age += dt_years * cell.pace

        cell.somatic_mutations += _poisson(
            self.rng, p.mutations_per_year.value * dt_years * self.env.mutagen_factor
        )

        if not cell.senescent:
            n_div = _poisson(self.rng, p.divisions_per_year.value * dt_years * self.env.proliferation_factor)
            keep = 1.0 - p.telomerase_compensation.value
            for _ in range(n_div):
                loss = max(
                    0.0,
                    self.rng.gauss(
                        p.telomere_loss_per_division_bp.value, p.telomere_loss_per_division_bp.value * 0.3
                    ),
                )
                cell.telomere_bp -= loss * keep
            cell.divisions += n_div

            hazard = self.params["senescence_hazard_base"] * math.exp(
                (self.params["telomere_senescence_bp"] - cell.telomere_bp)
                / self.params["telomere_hazard_scale_bp"]
            )
            hazard += (
                self.params["senescence_hazard_base"]
                * (cell.somatic_mutations / self.params["mutation_senescence_threshold"]) ** 2
            )
            if self.rng.random() < 1.0 - math.exp(-hazard * dt_years):
                cell.senescent = True

    def simulate_tissue(
        self,
        cell_type: str,
        years: float,
        n_cells: int = 200,
        dt_years: float = 0.5,
        report_every_years: float = 10.0,
    ) -> list[TissueReport]:
        cells = [self.new_cell(cell_type) for _ in range(n_cells)]
        reports: list[TissueReport] = []
        steps = int(round(years / dt_years))
        every = max(1, int(round(report_every_years / dt_years)))
        reports.append(self._report(cells, 0.0))
        for step in range(1, steps + 1):
            for c in cells:
                self.tick(c, dt_years)
            if step % every == 0 or step == steps:
                reports.append(self._report(cells, step * dt_years))
        return reports

    @staticmethod
    def _report(cells: list[CellState], age: float) -> TissueReport:
        n = len(cells)
        return TissueReport(
            age_years=age,
            n_cells=n,
            senescent_fraction=sum(c.senescent for c in cells) / n,
            mean_telomere_bp=sum(c.telomere_bp for c in cells) / n,
            mean_mutations=sum(c.somatic_mutations for c in cells) / n,
            mean_epigenetic_age=sum(c.epigenetic_age for c in cells) / n,
        )

    @staticmethod
    def evidence_table(cell_type: str) -> list[Parameter]:
        p = CELL_TYPES[cell_type]
        return [
            p.divisions_per_year,
            p.telomere_loss_per_division_bp,
            p.mutations_per_year,
            p.telomerase_compensation,
            *GLOBAL_PARAMS.values(),
        ]
