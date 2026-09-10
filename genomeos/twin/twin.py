"""BioTwin (task 3.3): one genome + one measured state + one environment.

A Twin is a JSON document under data/twins/. `fork` copies it with a change;
`run` advances its cell populations with the ageing runtime starting from the
*measured* state (telomere length, epigenetic age) rather than from birth;
`diff` compares two runs clock by clock. Variant consequences on timer genes
are translated into parameter modifiers with `inferred` evidence.
"""

from __future__ import annotations

import copy
import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

from genomeos.ir import Evidence, EvidenceKind, Parameter
from genomeos.lib import LIBRARIES
from genomeos.runtime.cell import CELL_TYPES, CellRuntime, CellTypeParams, Environment, TissueReport
from genomeos.runtime.uncertainty import UncertaintyReport, report_for_ageing

E_TWIN = Evidence(EvidenceKind.INFERRED, "GenomeOS twin rule: loss-of-function in a timer library gene")

# consequence -> how a loss-of-function in a library's gene changes the ageing parameters
LOF = {"nonsense", "frameshift_variant", "start_lost", "splice_site"}
TIMER_MODIFIERS: dict[str, dict[str, float]] = {
    "timer.telomere": {"telomerase_compensation": 0.0, "telomere_loss_per_division_bp": 1.2},
    "core.dna_repair": {"mutations_per_year": 1.5},
    "timer.progeroid_maintenance": {"mutations_per_year": 1.5, "divisions_per_year": 1.2},
    "timer.senescence_checkpoint": {"senescence_hazard_base": 0.5},
    "timer.growth_longevity_axis": {"epigenetic_pace": 0.9},
}


@dataclass(slots=True)
class MeasuredState:
    chronological_age: float = 0.0
    telomere_bp: float | None = None  # from reads (genomeos telomere) or assay
    epigenetic_age: float | None = None  # from methylation (genomeos clock)
    notes: str = ""


@dataclass(slots=True)
class Variant:
    gene: str
    consequence: str
    hgvs: str = ""
    source: str = ""


@dataclass(slots=True)
class Twin:
    name: str
    genome: str = ""  # reference FASTA or haplotype FASTA path
    vcf: str = ""
    sex: str = "unknown"
    measured: MeasuredState = field(default_factory=MeasuredState)
    environment: Environment = field(default_factory=Environment)
    variants: list[Variant] = field(default_factory=list)
    modifiers: dict[str, float] = field(
        default_factory=dict
    )  # parameter -> multiplier (or absolute for compensation)
    parent: str = ""
    note: str = ""

    # ---- persistence -------------------------------------------------------

    def save(self, directory: str | Path = "data/twins") -> Path:
        d = Path(directory)
        d.mkdir(parents=True, exist_ok=True)
        p = d / f"{self.name}.json"
        p.write_text(json.dumps(asdict(self), indent=2))
        return p

    @classmethod
    def load(cls, path: str | Path) -> Twin:
        d = json.loads(Path(path).read_text())
        d["measured"] = MeasuredState(**d.get("measured", {}))
        d["environment"] = Environment(**d.get("environment", {}))
        d["variants"] = [Variant(**v) for v in d.get("variants", [])]
        return cls(**d)

    def fork(self, name: str, note: str = "") -> Twin:
        t = copy.deepcopy(self)
        t.name, t.parent, t.note = name, self.name, note
        return t

    # ---- variants -> parameters -------------------------------------------

    def add_variant(self, gene: str, consequence: str, hgvs: str = "", source: str = "") -> list[str]:
        self.variants.append(Variant(gene, consequence, hgvs, source))
        applied = timer_effects_for_variants([self.variants[-1]])
        for k, v in applied.items():
            self.modifiers[k] = v if k in ("telomerase_compensation",) else self.modifiers.get(k, 1.0) * v
        return sorted(applied)

    # ---- run -------------------------------------------------------------

    def run(
        self,
        cell_type: str = "fibroblast",
        years: float = 40.0,
        cells: int = 500,
        seed: int = 1,
        use_calibration: bool = True,
    ) -> TwinRun:
        base = CELL_TYPES[cell_type]
        if use_calibration:
            base = _calibrated(base)
        params = _modified_params(base, self.modifiers)
        env = copy.copy(self.environment)
        if "epigenetic_pace" in self.modifiers:
            env.epigenetic_pace *= self.modifiers["epigenetic_pace"]
        rt = CellRuntime(seed=seed, env=env)
        if "senescence_hazard_base" in self.modifiers:
            rt.params["senescence_hazard_base"] *= self.modifiers["senescence_hazard_base"]
        rt_types = dict(CELL_TYPES)
        rt_types[cell_type] = params
        cellsl = []
        for _ in range(cells):
            c = rt.new_cell(cell_type)
            c.age_years = self.measured.chronological_age
            if self.measured.telomere_bp is not None:
                c.telomere_bp = max(1000.0, rt.rng.gauss(self.measured.telomere_bp, 400.0))
            c.epigenetic_age = (
                self.measured.epigenetic_age
                if self.measured.epigenetic_age is not None
                else self.measured.chronological_age
            )
            c.somatic_mutations = int(base.mutations_per_year.value * self.measured.chronological_age)
            cellsl.append(c)
        reports = [rt._report(cellsl, self.measured.chronological_age)]
        dt = 0.5
        steps = int(round(years / dt))
        import genomeos.runtime.cell as cellmod

        original = cellmod.CELL_TYPES[cell_type]
        cellmod.CELL_TYPES[cell_type] = params
        try:
            for step in range(1, steps + 1):
                for c in cellsl:
                    rt.tick(c, dt)
                if step % 2 == 0 or step == steps:
                    reports.append(rt._report(cellsl, self.measured.chronological_age + step * dt))
        finally:
            cellmod.CELL_TYPES[cell_type] = original
        evidence = [
            params.divisions_per_year,
            params.telomere_loss_per_division_bp,
            params.mutations_per_year,
            params.telomerase_compensation,
            *rt.evidence_table(cell_type)[4:],
        ]
        return TwinRun(self.name, cell_type, reports, report_for_ageing(evidence), dict(self.modifiers))


def _calibrated(base: CellTypeParams) -> CellTypeParams:
    """Apply a saved attrition calibration (genomeos calibrate) for this cell type, if any."""
    from genomeos.results import load_result

    r = load_result(f"calibration_{base.name}_attrition")
    if not r:
        return base
    ev = r.get("evidence", {})
    fitted = Parameter(
        "divisions_per_year",
        float(r["fitted_divisions_per_year"]),
        "1/yr",
        Evidence(EvidenceKind.INFERRED, ev.get("source", "calibration"), note=ev.get("note", "")),
        float(r.get("confidence", 0.4)),
    )
    return CellTypeParams(
        base.name,
        fitted,
        base.telomere_loss_per_division_bp,
        base.mutations_per_year,
        base.telomerase_compensation,
    )


def _modified_params(base: CellTypeParams, modifiers: dict[str, float]) -> CellTypeParams:
    def mod(p: Parameter, key: str, absolute: bool = False) -> Parameter:
        if key not in modifiers:
            return p
        value = modifiers[key] if absolute else p.value * modifiers[key]
        return Parameter(p.name, value, p.unit, E_TWIN, min(p.confidence, 0.4))

    return CellTypeParams(
        base.name,
        mod(base.divisions_per_year, "divisions_per_year"),
        mod(base.telomere_loss_per_division_bp, "telomere_loss_per_division_bp"),
        mod(base.mutations_per_year, "mutations_per_year"),
        mod(base.telomerase_compensation, "telomerase_compensation", absolute=True),
    )


def timer_effects_for_variants(variants: list[Variant]) -> dict[str, float]:
    """Parameter modifiers implied by loss-of-function variants in timer/maintenance library genes."""
    out: dict[str, float] = {}
    for v in variants:
        if v.consequence not in LOF:
            continue
        for lib_id, mods in TIMER_MODIFIERS.items():
            if v.gene in LIBRARIES[lib_id].genes:
                for k, val in mods.items():
                    out[k] = val if k == "telomerase_compensation" else out.get(k, 1.0) * val
    return out


@dataclass(slots=True)
class TwinRun:
    twin: str
    cell_type: str
    reports: list[TissueReport]
    uncertainty: UncertaintyReport
    modifiers: dict[str, float]

    @property
    def final(self) -> TissueReport:
        return self.reports[-1]

    def to_dict(self) -> dict:
        return {
            "twin": self.twin,
            "cell_type": self.cell_type,
            "modifiers": self.modifiers,
            "reports": [asdict(r) for r in self.reports],
            "uncertainty": self.uncertainty.to_dict(),
        }


def diff_runs(a: TwinRun, b: TwinRun) -> dict[str, dict[str, float]]:
    """Clock-by-clock difference at the end of two runs (b minus a)."""
    fa, fb = a.final, b.final
    return {
        "telomere_bp": {
            "a": fa.mean_telomere_bp,
            "b": fb.mean_telomere_bp,
            "delta": fb.mean_telomere_bp - fa.mean_telomere_bp,
        },
        "senescent_fraction": {
            "a": fa.senescent_fraction,
            "b": fb.senescent_fraction,
            "delta": fb.senescent_fraction - fa.senescent_fraction,
        },
        "mutations": {
            "a": fa.mean_mutations,
            "b": fb.mean_mutations,
            "delta": fb.mean_mutations - fa.mean_mutations,
        },
        "epigenetic_age": {
            "a": fa.mean_epigenetic_age,
            "b": fb.mean_epigenetic_age,
            "delta": fb.mean_epigenetic_age - fa.mean_epigenetic_age,
        },
    }
