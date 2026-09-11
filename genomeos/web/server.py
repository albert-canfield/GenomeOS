"""GenomeOS light web UI: a stdlib HTTP server exposing the runtime as JSON.

    genomeos serve            # http://127.0.0.1:8765

Binds to localhost only. File access is restricted to the project's data/
directory. No external dependencies; the page is a single static HTML file.
"""

from __future__ import annotations

import json
import threading
import webbrowser
from functools import lru_cache
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

from genomeos import __version__
from genomeos.genome import Locus, Sequence, read_fasta
from genomeos.lang import BioLangError, parse
from genomeos.lib import LIBRARIES
from genomeos.runtime import (
    CELL_TYPES,
    STANDARD_CODE,
    VERTEBRATE_MITOCHONDRIAL_CODE,
    CellRuntime,
    Environment,
    NetworkRuntime,
    find_orfs,
)

STATIC = Path(__file__).parent / "static"
TABLES = {"standard": STANDARD_CODE, "mito": VERTEBRATE_MITOCHONDRIAL_CODE}


@lru_cache(maxsize=8)
def _read_fasta_cached(path: str) -> dict[str, Sequence]:
    return read_fasta(path)


class ApiError(Exception):
    def __init__(self, message: str, status: int = 400) -> None:
        super().__init__(message)
        self.status = status


class Api:
    """All endpoints as plain methods so they can be unit-tested without HTTP."""

    def __init__(self, root: Path) -> None:
        self.root = root.resolve()
        self.data_dir = self.root / "data"

    # ---- helpers -------------------------------------------------------

    def _safe(self, rel: str) -> Path:
        p = (self.root / rel).resolve()
        if self.data_dir not in p.parents and p != self.data_dir:
            raise ApiError("path outside data/ is not allowed", 403)
        if not p.is_file():
            raise ApiError(f"no such file: {rel}", 404)
        return p

    @staticmethod
    def _load(path: str) -> dict[str, Sequence]:
        return _read_fasta_cached(path)

    # ---- files -----------------------------------------------------------

    def files(self) -> dict:
        genomes, modules, vcfs = [], [], []
        if self.data_dir.is_dir():
            for p in sorted(self.data_dir.rglob("*")):
                if not p.is_file():
                    continue
                rel = str(p.relative_to(self.root))
                if p.name.endswith((".fa", ".fa.gz", ".fasta", ".fasta.gz", ".fna", ".fna.gz")):
                    genomes.append({"path": rel, "size": p.stat().st_size})
                elif p.suffix == ".bio":
                    modules.append({"path": rel, "size": p.stat().st_size})
                elif p.suffix in (".vcf",) or p.name.endswith(".vcf.gz"):
                    vcfs.append({"path": rel, "size": p.stat().st_size})
        return {"genomes": genomes, "modules": modules, "vcfs": vcfs, "version": __version__}

    def module_source(self, rel: str) -> dict:
        return {"path": rel, "source": self._safe(rel).read_text()}

    # ---- genome ----------------------------------------------------------

    def genome_info(self, rel: str) -> dict:
        path = self._safe(rel)
        rows = []
        for name, seq in self._load(str(path)).items():
            fwd, rev = seq.telomeric_repeats()
            rows.append(
                {
                    "chromosome": name,
                    "length_bp": len(seq),
                    "gc": round(seq.gc_content(), 4),
                    "n_fraction": round(seq.n_fraction(), 4),
                    "ttaggg_fwd": fwd,
                    "ttaggg_rev": rev,
                }
            )
        return {"path": rel, "sequences": rows, "total_bp": sum(r["length_bp"] for r in rows)}

    def genome_fetch(self, rel: str, locus: str) -> dict:
        path = self._safe(rel)
        seqs = self._load(str(path))
        try:
            loc = Locus.parse(locus)
        except Exception as e:  # noqa: BLE001
            raise ApiError(f"bad locus {locus!r}: {e}") from None
        if loc.chrom not in seqs:
            raise ApiError(f"no chromosome {loc.chrom!r} in {rel}", 404)
        if loc.length > 20000:
            raise ApiError("window limited to 20,000 bp")
        seq = seqs[loc.chrom][loc.start : loc.end]
        if loc.strand.value == "-":
            seq = seq.reverse_complement()
        return {"locus": str(loc), "sequence": str(seq), "gc": round(seq.gc_content(), 4)}

    def genome_orfs(self, rel: str, min_aa: int = 50, table: str = "standard", limit: int = 500) -> dict:
        path = self._safe(rel)
        if table not in TABLES:
            raise ApiError(f"table must be one of {sorted(TABLES)}")
        out = []
        for name, seq in self._load(str(path)).items():
            for orf in find_orfs(seq, chrom=name, min_aa=min_aa, table=TABLES[table]):
                out.append(
                    {
                        "locus": str(orf.locus),
                        "chrom": orf.locus.chrom,
                        "start": orf.locus.start,
                        "end": orf.locus.end,
                        "strand": orf.locus.strand.value,
                        "length_aa": orf.length_aa,
                        "protein": orf.protein,
                    }
                )
        out.sort(key=lambda o: -o["length_aa"])
        return {"count": len(out), "orfs": out[:limit], "truncated": len(out) > limit}

    # ---- program ---------------------------------------------------------

    def compile(self, source: str) -> dict:
        try:
            module = parse(source, name_hint="editor")
        except BioLangError as e:
            return {"ok": False, "error": str(e)}
        return {
            "ok": True,
            "module": module.to_dict(),
            "report": {
                "entities": len(module.entities),
                "rules": len(module.rules),
                "parameters": len(module.parameters),
                "confidence": module.confidence_report(),
                "unknown": [
                    {"id": u.id, "locus": str(u.locus) if u.locus else None} for u in module.unknowns()
                ],
                "weak_rules": [
                    {
                        "id": r.id,
                        "confidence": r.confidence,
                        "evidence": r.evidence.kind.value,
                        "source": r.evidence.source,
                    }
                    for r in module.rules
                    if r.confidence < 0.5
                ],
            },
        }

    def run(
        self,
        source: str,
        hours: float = 50.0,
        dt: float = 0.01,
        initial: dict | None = None,
        context: dict | None = None,
        seed: int | None = None,
        points: int = 600,
    ) -> dict:
        try:
            module = parse(source, name_hint="editor")
        except BioLangError as e:
            return {"ok": False, "error": str(e)}
        hours = min(max(float(hours), 0.1), 2000.0)
        dt = min(max(float(dt), 1e-4), 1.0)
        vm = NetworkRuntime(module, context=context or {}, seed=seed)
        steps = int(round(hours / dt))
        record_every = max(1, steps // points)
        traj = vm.run(
            hours=hours,
            dt=dt,
            initial={k: float(v) for k, v in (initial or {}).items()},
            record_every=record_every,
        )
        return {
            "ok": True,
            "module": module.name,
            "active_rules": len(vm.active_rules),
            "total_rules": len(module.rules),
            "times": traj.times,
            "species": traj.species,
            "levels": traj.levels,
            "peaks": {s: traj.peaks(s) for s in traj.species},
        }

    # ---- ageing ----------------------------------------------------------

    def age(
        self,
        cell_type: str = "fibroblast",
        years: float = 90,
        cells: int = 500,
        seed: int = 1,
        proliferation: float = 1.0,
        mutagen: float = 1.0,
        pace: float = 1.0,
    ) -> dict:
        if cell_type not in CELL_TYPES:
            raise ApiError(f"cell_type must be one of {sorted(CELL_TYPES)}")
        years = min(max(float(years), 1), 120)
        cells = min(max(int(cells), 10), 5000)
        env = Environment(
            proliferation_factor=float(proliferation),
            mutagen_factor=float(mutagen),
            epigenetic_pace=float(pace),
        )
        rt = CellRuntime(seed=int(seed), env=env)
        reports = rt.simulate_tissue(cell_type, years, n_cells=cells, dt_years=0.5, report_every_years=1.0)
        return {
            "cell_type": cell_type,
            "cells": cells,
            "years": years,
            "reports": [
                {
                    "age": r.age_years,
                    "senescent": r.senescent_fraction,
                    "telomere_bp": r.mean_telomere_bp,
                    "mutations": r.mean_mutations,
                    "epigenetic_age": r.mean_epigenetic_age,
                }
                for r in reports
            ],
            "evidence": [
                {
                    "name": p.name,
                    "value": p.value,
                    "unit": p.unit,
                    "confidence": p.confidence,
                    "kind": p.evidence.kind.value,
                    "source": p.evidence.source,
                    "note": p.evidence.note,
                }
                for p in rt.evidence_table(cell_type)
            ],
            "cell_types": sorted(CELL_TYPES),
            "uncertainty": __import__("genomeos.runtime.uncertainty", fromlist=["report_for_ageing"])
            .report_for_ageing(rt.evidence_table(cell_type))
            .to_dict(),
        }

    # ---- twins ------------------------------------------------------------

    def individuals(self) -> dict:
        """The test human and every imported person (local files only)."""
        from genomeos.genome.individuals import checks, list_individuals

        root = self.root / "data" / "individuals"
        people = list_individuals(root)
        for p in people:
            p["checks"] = checks(p["name"], root)
        return {"individuals": people}

    def individual_screen(self, name: str) -> dict:
        """ClinVar pathogenic alleles the person carries (ClinVar distilled once, locally)."""
        from genomeos.genome import clinvar

        knowledge = self.root / "data" / "knowledge" / "clinvar"
        if not clinvar.pathogenic_path(knowledge).exists():
            clinvar.distil(knowledge)
        try:
            return clinvar.screen(name, None, self.root / "data" / "individuals", knowledge)
        except FileNotFoundError as ex:
            raise ApiError(str(ex)) from ex

    def individual_check(self, name: str, chrom: str) -> dict:
        """The person's variants of one chromosome applied to the local reference: statistics and verdict."""
        from genomeos.genome.individuals import check

        try:
            return check(name, chrom, self.root / "data" / "individuals", self.root / "data" / "reference")
        except FileNotFoundError as ex:
            raise ApiError(str(ex)) from ex

    def individual_import(self, path: str, name: str, note: str = "", replace: bool = False) -> dict:
        """Split a VCF on this machine into per-chromosome files under data/individuals/<name>."""
        from genomeos.genome.individuals import import_vcf

        try:
            root = self.root / "data" / "individuals"
            return {"ok": True, "manifest": import_vcf(path, name, note, root, replace)}
        except (ValueError, FileNotFoundError, FileExistsError) as ex:
            raise ApiError(str(ex)) from ex

    def individual_genes(self, name: str, chrom: str, limit: int = 40) -> dict:
        """One chromosome gene by gene for one person: variants and coding consequences."""
        from genomeos.genome import Annotation, IndexedGenome, default_gencode
        from genomeos.genome.individuals import gene_by_gene

        gff = default_gencode({chrom})
        fa = self.root / "data" / "reference" / f"{chrom}.fa"
        if not gff or not fa.exists():
            raise ApiError(f"{chrom} needs local models and sequence (genomeos data fetch --chrom {chrom})")
        ann = Annotation.from_gff3(gff, {chrom})
        genome = IndexedGenome(str(fa))
        try:
            return gene_by_gene(name, chrom, ann, genome, self.root / "data" / "individuals", limit)
        except FileNotFoundError as ex:
            raise ApiError(str(ex)) from ex
        finally:
            genome.close()

    def individual_predict(self, name: str, gene: str, chrom: str, max_variants: int = 30) -> dict:
        """Feature d: the person's regulatory variants on one gene, predicted per variant (cached)."""
        from genomeos.genome import Annotation, IndexedGenome, default_gencode
        from genomeos.genome.individuals import vcf_path
        from genomeos.genome.regulation import regulation_of
        from genomeos.genome.regulatory import load_ccres
        from genomeos.predict import AlphaGenomeAdapter, status
        from genomeos.predict.individual_effects import predict_gene

        st = status()
        if not st["enabled"]:
            raise ApiError(f"AlphaGenome predictions are disabled: {st['reason']}")
        vcf = vcf_path(name, chrom, self.root / "data" / "individuals")
        if vcf is None:
            raise ApiError(f"{name} has no rows on {chrom}")
        gff = default_gencode({chrom})
        fa = self.root / "data" / "reference" / f"{chrom}.fa"
        ccres = load_ccres(chrom)
        if not gff or not fa.exists() or not ccres:
            raise ApiError(
                f"{chrom} needs local models, sequence and elements (genomeos data fetch --chrom {chrom})"
            )
        ann = Annotation.from_gff3(gff, {chrom})
        genome = IndexedGenome(str(fa))
        try:
            reg = regulation_of(gene.upper(), chrom, ccres, ann, genome.lengths[chrom])
        except KeyError as ex:
            raise ApiError(f"{gene} not on {chrom}") from ex
        finally:
            genome.close()
        scorer = AlphaGenomeAdapter()._live_scorer(threshold=0.02)  # noqa: SLF001
        r = predict_gene(name, reg["gene"], chrom, reg, scorer, vcf, max_variants=max_variants)
        return {"model": st["model"], **r}

    def twins(self) -> dict:
        from genomeos.twin import Twin

        d = self.root / "data" / "twins"
        out = []
        if d.is_dir():
            for p in sorted(d.glob("*.json")):
                try:
                    t = Twin.load(p)
                except Exception:  # noqa: BLE001
                    continue
                out.append(
                    {
                        "name": t.name,
                        "parent": t.parent,
                        "age": t.measured.chronological_age,
                        "telomere_bp": t.measured.telomere_bp,
                        "epigenetic_age": t.measured.epigenetic_age,
                        "variants": [f"{v.gene}:{v.consequence}" for v in t.variants],
                        "modifiers": t.modifiers,
                        "note": t.note,
                    }
                )
        return {"twins": out}

    def twin_fork(
        self,
        name: str,
        new_name: str,
        variants: list[str],
        mutagen: float | None,
        pace: float | None,
        note: str = "",
    ) -> dict:
        from genomeos.twin import Twin

        src = self.root / "data" / "twins" / f"{name}.json"
        if not src.is_file():
            raise ApiError(f"no twin {name!r}", 404)
        if not new_name.replace("_", "").replace("-", "").isalnum():
            raise ApiError("new name must be alphanumeric with _ or -")
        t = Twin.load(src).fork(new_name, note)
        applied = {}
        for spec in variants:
            gene, _, cons = spec.partition(":")
            if gene:
                applied[spec] = t.add_variant(gene.strip(), (cons or "nonsense").strip())
        if mutagen is not None:
            t.environment.mutagen_factor = float(mutagen)
        if pace is not None:
            t.environment.epigenetic_pace = float(pace)
        t.save(self.root / "data" / "twins")
        return {"ok": True, "twin": new_name, "applied": applied, "modifiers": t.modifiers}

    def twin_new(
        self,
        name: str,
        age: float,
        telomere: float | None,
        epigenetic_age: float | None,
        sex: str = "unknown",
    ) -> dict:
        from genomeos.twin import Twin
        from genomeos.twin.twin import MeasuredState

        if not name.replace("_", "").replace("-", "").isalnum():
            raise ApiError("name must be alphanumeric with _ or -")
        t = Twin(name, sex=sex, measured=MeasuredState(float(age), telomere, epigenetic_age))
        t.save(self.root / "data" / "twins")
        return {"ok": True, "twin": name}

    def twin_run(
        self,
        names: list[str],
        cell_type: str = "fibroblast",
        years: float = 40,
        cells: int = 300,
        seed: int = 1,
    ) -> dict:
        from genomeos.twin import Twin, diff_runs

        if cell_type not in CELL_TYPES:
            raise ApiError(f"cell_type must be one of {sorted(CELL_TYPES)}")
        runs = []
        for name in names[:2]:
            src = self.root / "data" / "twins" / f"{name}.json"
            if not src.is_file():
                raise ApiError(f"no twin {name!r}", 404)
            t = Twin.load(src)
            runs.append(
                t.run(cell_type, years=min(float(years), 100), cells=min(int(cells), 3000), seed=int(seed))
            )
        out = {"runs": [r.to_dict() for r in runs]}
        if len(runs) == 2:
            out["diff"] = diff_runs(runs[0], runs[1])
        return out

    # ---- organism (grow) -------------------------------------------------

    def grow(self, module: str, until: str = "800", depth: int = 3) -> dict:
        """Grow a BioLang v0.3 organism program and score it against its reference when one exists."""
        from genomeos.ir import to_minutes
        from genomeos.lang import parse_file
        from genomeos.organism import REFERENCES, ReferenceLineage
        from genomeos.organism.diff import compare
        from genomeos.runtime.body import Body

        path = self._safe(module)
        m = parse_file(path)
        if m.organism is None:
            raise ApiError(f"{module} has no organism block")
        parts = str(until).split()
        until_min = to_minutes(float(parts[0]), parts[1] if len(parts) > 1 else "min")
        until_min = min(until_min, to_minutes(100, "yr"))
        body = Body(m, max_cells=50_000).run(until=until_min)
        marks = [0.5, 25, 50, 100, 150, 200, 250, 300, 350, 400, 500, 600, 800, 1500, 2000, 3000, 4000, 5700]
        if until_min > 10_000:
            days = (1, 3, 5, 14, 21, 56, 100, 180, 266, 365, 730, 1826, 3652, 6574, 7305, 10957, 14610, 29220)
            marks = [to_minutes(d, "d") for d in days]
        marks = [t for t in marks if t <= until_min] + [until_min]
        ref = None
        if m.organism.reference in REFERENCES and Path(REFERENCES[m.organism.reference]).exists():
            ref = ReferenceLineage.load(REFERENCES[m.organism.reference])
        out = {
            "summary": body.summary(),
            "organism": {
                "name": m.organism.name,
                "species": m.organism.species,
                "resolution": m.organism.resolution,
                "factors": m.organism.factors,
                "decisions": len(m.decisions),
                "timers": len(m.timers),
                "signals": len(m.signals()),
                "stages": [
                    {"name": st.name, "start": st.start, "end": st.end, "unit": st.unit} for st in m.stages
                ],
            },
            "curve": [
                {
                    "t": t,
                    "cells": body.count_at(t),
                    "deaths": body.deaths_by(t),
                    "reference": ref.count_at(t) if ref else None,
                }
                for t in marks
            ],
            "fates": body.fates_at(until_min),
            "tree": body.tree(depth=max(0, min(int(depth), 6))),
            "space": {
                "width": m.organism.width,
                "height": m.organism.height,
                "rows": body.type_map(),
                "bands": body.bands_along_x(),
                "blocked": body.blocked_divisions,
            }
            if body.spatial
            else None,
            "asserts": body.check_asserts(),
            "uncertainty": body.uncertainty().to_dict(),
            "unknown": dict(body.unknown),
            "fired": [
                {"id": k, "n": v}
                for k, v in body.fired.most_common()
                if not k.startswith(("div_", "fate_", "die_"))
            ][:40],
        }
        if ref is not None:
            out["diff"] = compare(body, ref, until=until_min).to_dict()
        return out

    # ---- development (space) --------------------------------------------

    def develop(self, model: str, width: int = 60, hours: float = 100.0) -> dict:
        width = min(max(int(width), 20), 200)
        hours = min(max(float(hours), 5), 400)
        if model == "flag":
            from genomeos.runtime.spatial import french_flag

            rt = french_flag(width=width, height=6, hours=hours)
            return {
                "model": "flag",
                "rows": rt.type_map(),
                "bands": rt.bands_along_x(),
                "census": rt.census(),
                "profile": rt.fields["morphogen"].profile_x(),
                "parameters": [
                    {
                        "name": q.name,
                        "value": q.value,
                        "kind": q.evidence.kind.value,
                        "source": q.evidence.source,
                        "confidence": q.confidence,
                    }
                    for q in rt.parameters
                ],
            }
        if model == "segmentation":
            from genomeos.runtime.segmentation import run_segmentation

            r = run_segmentation(length=width, hours=hours)
            cells = ["A" if ph < 3.14159265 else "P" for ph in r.frozen_phase if ph == ph]
            return {
                "model": "segmentation",
                "period_h": r.period_h,
                "count": r.count,
                "expected": r.expected,
                "segments": r.segments,
                "cells": cells,
                "fgf": r.fields["fgf"],
                "uncertainty": r.uncertainty().to_dict(),
            }
        if model == "gastrulation":
            from genomeos.runtime.gastrulation import run_gastrulation

            r = run_gastrulation(cells=min(width, 120), hours=min(hours, 60))
            return {
                "model": "gastrulation",
                "fates": r.fates,
                "positions": r.positions,
                "proportions": r.proportions(),
                "uncertainty": r.uncertainty().to_dict(),
            }
        raise ApiError("model must be flag, segmentation or gastrulation")

    # ---- debugger -----------------------------------------------------------

    def debug(
        self,
        source: str,
        breakpoints: list[str],
        initial: dict | None,
        until: float = 100.0,
        explain: list[str] | None = None,
        cell_type: str | None = None,
    ) -> dict:
        from genomeos.runtime.debugger import AgeingDebugger, NetworkDebugger

        if cell_type:
            if cell_type not in CELL_TYPES:
                raise ApiError(f"cell_type must be one of {sorted(CELL_TYPES)}")
            dbg = AgeingDebugger(cell_type, seed=1)
            for b in breakpoints:
                try:
                    dbg.add_breakpoint(b)
                except ValueError as e:
                    raise ApiError(str(e)) from None
            hit = dbg.step(years=min(float(until), 120))
            return {
                "ok": True,
                "kind": "ageing",
                "hit": str(hit) if hit else None,
                "time": dbg.session.time,
                "state": dbg.session.state,
                "trace": [
                    {
                        "time": line.time,
                        "subject": line.subject,
                        "message": line.message,
                        "evidence": line.evidence,
                        "confidence": line.confidence,
                    }
                    for line in dbg.session.trace[-40:]
                ],
                "explain": [
                    {
                        "time": line.time,
                        "subject": line.subject,
                        "message": line.message,
                        "evidence": line.evidence,
                        "confidence": line.confidence,
                    }
                    for line in dbg.explain()
                ],
            }
        try:
            module = parse(source, name_hint="editor")
        except BioLangError as e:
            return {"ok": False, "error": str(e)}
        dbg = NetworkDebugger(module)
        if initial:
            dbg.set_initial({k: float(v) for k, v in initial.items()})
        for b in breakpoints:
            try:
                dbg.add_breakpoint(b)
            except ValueError as e:
                raise ApiError(str(e)) from None
        hit = dbg.run_until_break(max_hours=min(float(until), 1000))
        targets = explain or ([hit.variable] if hit else [])
        expl = {
            sp: [
                {
                    "time": line.time,
                    "subject": line.subject,
                    "message": line.message,
                    "evidence": line.evidence,
                    "confidence": line.confidence,
                }
                for line in dbg.explain(sp)
            ]
            for sp in targets
            if sp in dbg.session.state
        }
        return {
            "ok": True,
            "kind": "network",
            "hit": str(hit) if hit else None,
            "time": dbg.session.time,
            "state": dbg.session.state,
            "explain": expl,
        }

    # ---- measurements (distilled results) -------------------------------

    def results(self) -> dict:
        from genomeos.results import list_results, load_result

        out = []
        for r in list_results(self.root / "data" / "results"):
            d = load_result(r["name"], self.root / "data" / "results") or {}
            out.append(
                {
                    "name": r["name"],
                    "date": r["date"],
                    "size": r["size"],
                    "summary": {
                        k: v
                        for k, v in d.items()
                        if isinstance(v, int | float | str) and k not in ("result",)
                    },
                }
            )
        return {"results": out}

    def twin_measure(self, name: str) -> dict:
        """Apply the streamed HG002 measurements (telomere, methylation) to a twin's measured state."""
        from genomeos.results import load_result
        from genomeos.twin import Twin

        src = self.root / "data" / "twins" / f"{name}.json"
        if not src.is_file():
            raise ApiError(f"no twin {name!r}", 404)
        t = Twin.load(src)
        applied = {}
        tel = load_result("hg002_telomere_stream", self.root / "data" / "results")
        if tel:
            t.measured.telomere_bp = float(tel["southern_equivalent_bp_inferred"])
            applied["telomere_bp"] = t.measured.telomere_bp
        meth = load_result("hg002_methylation_stream", self.root / "data" / "results")
        if meth:
            t.measured.epigenetic_age = float(meth["horvath_age"])
            applied["epigenetic_age"] = t.measured.epigenetic_age
        t.measured.notes = (
            "measured state from streamed GIAB data: " + ", ".join(applied) if applied else t.measured.notes
        )
        t.save(self.root / "data" / "twins")
        return {"ok": True, "twin": name, "applied": applied}

    # ---- anatomy ----------------------------------------------------------

    def anatomy(self, rel: str | None, chrom: str | None) -> dict:
        from genomeos.genome import Annotation, anatomy_of, default_gencode, design_lessons
        from genomeos.results import load_result

        if not rel:
            r = load_result("anatomy_comparison", self.root / "data" / "results")
            if not r:
                raise ApiError("no saved comparison", 404)
            return r
        path = self._safe(rel)
        seqs = self._load(str(path))
        chrom = chrom or next(iter(seqs))
        if chrom not in seqs:
            raise ApiError(f"no chromosome {chrom!r}", 404)
        ann = None
        gff = default_gencode({chrom}) if chrom.startswith("chr") else None
        if "celegans" in rel:
            g = self.root / "data" / "reference" / "celegans" / "WBcel235.63.gff3.gz"
            gff = g if g.exists() else None
        if gff:
            ann = Annotation.from_gff3(gff, {chrom})
        a = anatomy_of(f"{path.name} {chrom}", seqs[chrom], ann, chrom if ann else None)
        d = a.to_dict()
        d["lessons"] = design_lessons([a]) if ann else []
        return d

    # ---- progress: jobs, project log ------------------------------------

    def jobs(self) -> dict:
        from genomeos import jobs

        return {"jobs": [j.to_dict() for j in jobs.all_status(self.root)]}

    def job_start(self, name: str) -> dict:
        from genomeos import jobs

        try:
            return {"ok": True, "job": jobs.start(name, self.root).to_dict()}
        except KeyError as e:
            raise ApiError(str(e)) from None

    def progress(self) -> dict:
        """The project's progress table (docs/PROGRESS.md) as rows."""
        p = self.root / "docs" / "PROGRESS.md"
        rows = []
        if p.exists():
            for line in p.read_text().splitlines():
                if (
                    not line.startswith("|")
                    or set(line.strip("| ")) <= {"-", " ", "|"}
                    or line.startswith("| Task")
                ):
                    continue
                cells = [c.strip() for c in line.strip().strip("|").split("|")]
                if len(cells) >= 3:
                    rows.append({"task": cells[0], "status": cells[1], "evidence": cells[2]})
        return {
            "rows": rows,
            "count": len(rows),
            "done": sum(1 for r in rows if r["status"].startswith("done")),
        }

    def unknown_wide(self) -> dict:
        """Genome-wide UNKNOWN classification: per-chromosome rows and class totals, as far as it has run."""
        from genomeos.results import load_result

        r = load_result("unknown_genome_wide", self.root / "data" / "results") or {"chromosomes": {}}
        ch = r.get("chromosomes", {})
        order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
        rows = []
        for c in order:
            v = ch.get(c) or (
                load_result(f"unknown_{c}", self.root / "data" / "results") if c == "chr21" else None
            )
            if not v:
                continue
            by = {k: (x["bp"] if isinstance(x, dict) else x) for k, x in v["by_class"].items()}
            top = max(
                (kv for kv in by.items() if kv[0] != "unclassified"), key=lambda kv: kv[1], default=("-", 0)
            )
            if "curated_repeats" not in v:
                per = load_result(f"unknown_{c}", self.root / "data" / "results") or {}
                v = {**v, "curated_repeats": per.get("curated_repeats", False)}
            rows.append(
                {
                    "chrom": c,
                    "blocks": v["unknown_blocks"],
                    "unknown_mb": round(v["unknown_bp"] / 1e6, 1),
                    "classified": v["classified_fraction"],
                    "top_class": top[0],
                    "top_share": round(top[1] / v["unknown_bp"], 3) if v["unknown_bp"] else None,
                    "curated_repeats": bool(v.get("curated_repeats")),
                    "seconds": v.get("seconds"),
                }
            )
        tot_bp = sum(v["unknown_bp"] for v in ch.values()) or 1
        by_all: dict[str, int] = {}
        for v in ch.values():
            for k, x in v["by_class"].items():
                by_all[k] = by_all.get(k, 0) + (x["bp"] if isinstance(x, dict) else x)
        return {
            "done": len(rows),
            "total": 24,
            "unknown_bp": tot_bp,
            "classified": round(1 - by_all.get("unclassified", 0) / tot_bp, 4) if ch else None,
            "by_class": dict(sorted(by_all.items(), key=lambda kv: -kv[1])[:12]),
            "table": rows,
        }

    def individual_wide(self) -> dict:
        """The test human on every autosome (statistics only) and the reader on every chromosome."""
        from genomeos.results import load_result

        rd = self.root / "data" / "results"
        return {
            "twin": load_result("hg002_twin_by_chromosome", rd),
            "reader": load_result("reader_genome_wide", rd),
            "enhancer_targets": {
                p.stem.split("_")[-1]: load_result(p.stem, rd)
                for p in sorted(rd.glob("enhancer_targets_chr*.json"))
            },
            "enhancer_targets_genome_wide": load_result("enhancer_targets_genome_wide", rd),
            "mouse": {
                p.stem.split("_")[-1]: {
                    k: v
                    for k, v in (load_result(p.stem, rd) or {}).items()
                    if k not in ("node_rows", "mouse_domains")
                }
                for p in sorted(rd.glob("mouse_mm10_chr*.json"))
            },
        }

    def proteome_summary(self) -> dict:
        from genomeos.results import load_result

        out = {}
        for c in [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]:
            r = load_result(f"proteome_{c}", self.root / "data" / "results")
            if r:
                out[c] = {
                    "coding_genes": r["coding_genes"],
                    "coverage": r["coverage_fraction"],
                    "no_entry": len(r.get("no_reviewed_entry", [])),
                }
        return out

    def genome_wide(self) -> dict:
        from genomeos.results import load_result

        r = load_result("anatomy_hg38_by_chromosome", self.root / "data" / "results")
        if not r:
            raise ApiError("no genome-wide inventory yet", 404)
        chroms = r["chromosomes"]
        table = []
        for name, d in chroms.items():
            f = d["composition_fraction"]
            table.append(
                {
                    "chrom": name,
                    "length_mb": round(d["length"] / 1e6, 1),
                    "coding_genes": d["genes"]["protein_coding"],
                    "genes_per_mb": d["layout"]["coding_genes_per_mb"],
                    "cds": f["cds"],
                    "intron": f["intron"],
                    "intergenic": f["intergenic"],
                    "gap": f["gap"],
                    "exons": d["structure"]["exons_per_transcript_median"],
                    "intron_median": d["structure"]["intron_length_median"],
                    "cpg_islands": d["elements"]["cpg_islands"],
                    "seconds": d.get("seconds"),
                }
            )
        total_len = sum(d["length"] for d in chroms.values())
        totals = {
            k: sum(d["composition_bp"].get(k, 0) for d in chroms.values())
            for k in ("cds", "intron", "intergenic", "gap", "noncoding_exon_or_utr")
        }
        return {
            "assembly": r.get("assembly"),
            "done": len(chroms),
            "total": 25,
            "table": table,
            "genome_bp": total_len,
            "coding_genes": sum(d["genes"]["protein_coding"] for d in chroms.values()),
            "composition_bp": totals,
        }

    # ---- blocks (2-D viewer/editor) ---------------------------------------

    def _annotation_for(self, rel: str, chrom: str):
        from genomeos.genome import Annotation, default_gencode

        key = (rel, chrom)
        if not hasattr(self, "_ann_cache"):
            self._ann_cache = {}
        if key in self._ann_cache:
            return self._ann_cache[key]
        gff = default_gencode({chrom}) if chrom.startswith("chr") else None
        if "celegans" in rel:
            g = self.root / "data" / "reference" / "celegans" / "WBcel235.63.gff3.gz"
            gff = g if g.exists() else None
        ann = Annotation.from_gff3(gff, {chrom}) if gff else None
        if len(self._ann_cache) > 6:
            self._ann_cache.pop(next(iter(self._ann_cache)))
        self._ann_cache[key] = ann
        return ann

    def blocks(self, rel: str, chrom: str, start: int, end: int) -> dict:
        from genomeos.genome.blocks import ELEMENT_WINDOW, blocks_for_window

        path = self._safe(rel)
        seqs = self._load(str(path))
        chrom = chrom or next(iter(seqs))
        if chrom not in seqs:
            raise ApiError(f"no chromosome {chrom!r}", 404)
        length = len(seqs[chrom])
        start = max(0, int(start))
        end = min(length, int(end)) if end else length
        if end <= start:
            raise ApiError("empty window")
        ann = self._annotation_for(rel, chrom)
        if ann is None:
            from genomeos.genome import Annotation

            ann = Annotation()
        seq = seqs[chrom][start:end] if end - start <= ELEMENT_WINDOW else None
        return blocks_for_window(ann, chrom, start, end, seq, length)

    def blocks_check_move(
        self, blocks: list, block_id: str, new_start: int, chrom_length: int | None
    ) -> dict:
        from genomeos.genome.blocks import check_move

        return check_move(blocks, block_id, int(new_start), chrom_length)

    # ---- cancer -------------------------------------------------------------

    def cancer_knowledge(self, symbol: str | None = None, top: int = 30) -> dict:
        from genomeos.results import load_result

        k = load_result("cancer_msk_impact_2017", self.root / "data" / "results")
        if not k:
            raise ApiError("no distilled cancer knowledge; run genomeos cancer distil", 404)
        if symbol:
            v = k["genes"].get(symbol)
            if not v:
                raise ApiError(f"{symbol} not in the distilled panel", 404)
            return {"gene": symbol, "study": k["study"], "samples": k["samples"], **v}
        rows = sorted(k["genes"].items(), key=lambda kv: -kv[1]["frequency"])[:top]
        return {
            "study": k["study"],
            "samples": k["samples"],
            "cancer_types": k["cancer_types"],
            "genes": [
                {
                    "gene": g,
                    "frequency": v["frequency"],
                    "samples_mutated": v["samples_mutated"],
                    "hotspots": v["hotspots"][:3],
                    "top_types": sorted(v["by_cancer_type"].items(), key=lambda x: -x[1])[:4],
                }
                for g, v in rows
            ],
            "evidence": k["evidence"],
        }

    def cancer_tumour(self, vcf: str, deep: int = 5) -> dict:
        """Tumour-only pipeline on a local VCF; the packet is included."""
        from genomeos.cancer import analyse, tumour_packet

        path = self._safe(vcf)
        a = analyse(str(path), deep=deep)
        a["packet"] = tumour_packet(a)
        return a

    def therapeutics(
        self,
        vcf: str,
        hla: str = "",
        rna: str = "",
        normal: str = "",
        top: int = 8,
        offline: bool = False,
        cohort: str = "",
    ) -> dict:
        """Therapeutic target and mechanism reasoning for one tumour VCF."""
        from genomeos.therapeutics import analyse_vcf, machine_report, text_report
        from genomeos.therapeutics.design import dataset, design_readiness, negative_targets

        path = self._safe(vcf)
        a = analyse_vcf(
            str(path),
            normal_vcf=str(self._safe(normal)) if normal else None,
            hla=[x for x in hla.split(",") if x.strip()],
            rna=str(self._safe(rna)) if rna else None,
            top_genes=top,
            cohort=cohort,
            net=not offline,
        )
        out = machine_report(a)
        out["report"] = text_report(a)
        out["data_level_detail"] = a["data_level"]
        out["design"] = {
            c.gene: {
                "readiness": design_readiness(c),
                "negative_targets": negative_targets(c, a.get("provider_bundle")),
                "class_reason": c.class_reason,
                "localization": c.localization.to_dict(),
                "trafficking": c.trafficking.to_dict(),
                "normal_tissue": c.normal_tissue.to_dict(),
                "structure": c.structure.to_dict(),
                "neoantigen": c.neoantigen.to_dict() if c.neoantigen else None,
                "scores": c.scores.to_dict(),
                "mechanisms": [m.to_dict() for m in c.therapeutic_mechanisms],
            }
            for c in a["candidates"]
        }
        out["schema"] = dataset(a, None)["schema"]
        return out

    def cancer_compare(self, normal: str, tumour: str, genome: str, chrom: str) -> dict:
        from genomeos.cancer import agent_packet, annotate, somatic, suggest_cancer_type, surface_targets
        from genomeos.genome import Annotation, IndexedGenome, default_gencode
        from genomeos.results import load_result

        k = load_result("cancer_msk_impact_2017", self.root / "data" / "results") or {}
        n, t, g = self._safe(normal), self._safe(tumour), self._safe(genome)
        chroms = {chrom} if chrom else None
        som = somatic(str(n), str(t), chroms)
        gff = default_gencode(chroms)
        if not gff:
            raise ApiError("no annotation for that chromosome (chr21/chrM available)")
        ann = Annotation.from_gff3(gff, chroms)
        idx = IndexedGenome(g)
        ranked = annotate(som, ann, idx, k)
        idx.close()
        genes = {
            s.gene
            for s in ranked
            if s.gene and s.consequence not in ("synonymous_variant", "intron_variant", "intergenic")
        }
        types = suggest_cancer_type(genes, k) if k else []
        targets = surface_targets(genes)
        return {
            "somatic": len(som),
            "ranked": [s.to_dict() for s in ranked[:60]],
            "types": types,
            "targets": targets,
            "packet": agent_packet(str(t), ranked, types, targets),
        }

    # ---- molecules ------------------------------------------------------------

    def protein(self, gene: str, chrom: str | None) -> dict:
        from genomeos.genome import IndexedGenome, default_gencode
        from genomeos.molecules import protein_report

        if not gene or not gene.replace("-", "").replace("_", "").isalnum():
            raise ApiError("gene symbol required")
        ann = genome = None
        if chrom:
            gff = default_gencode({chrom})
            if gff:
                ann = self._annotation_for(f"data/reference/{chrom}.fa.gz", chrom)
                fa = self.root / "data" / "reference" / f"{chrom}.fa.gz"
                genome = IndexedGenome(fa) if fa.exists() else None
        try:
            return protein_report(gene.upper(), ann, genome, structure=True, coords=True)
        finally:
            if genome:
                genome.close()

    def report(self, gene: str, chrom: str) -> dict:
        """The gene dossier (Markdown and the dict behind it)."""
        from genomeos.genome import IndexedGenome, default_gencode
        from genomeos.report import gene_report, to_markdown

        if not gene or not gene.replace("-", "").replace("_", "").isalnum():
            raise ApiError("gene symbol required")
        fa = self.root / "data" / "reference" / f"{chrom}.fa.gz"
        if not default_gencode({chrom}) or not fa.exists():
            raise ApiError(f"{chrom} is not local; fetch it first")
        ann = self._annotation_for(f"data/reference/{chrom}.fa.gz", chrom)
        genome = IndexedGenome(fa)
        try:
            rep = gene_report(gene.upper(), chrom, ann, genome, self.root)
        except KeyError as e:
            raise ApiError(f"{gene} not on {chrom}") from e
        finally:
            genome.close()
        rep["markdown"] = to_markdown(rep)
        return rep

    def features(self) -> list[dict]:
        """Optional features (external keys or packages): always listed, disabled without them."""
        from genomeos.predict import status

        return [status()]

    def predict(self, variant: str) -> dict:
        """AlphaGenome predicted tissue effects for one variant, or the disabled status."""
        import re

        from genomeos.predict import AlphaGenomeAdapter, status

        st = status()
        if not st["enabled"]:
            return {"enabled": False, "reason": st["reason"], "how": st["how"]}
        m = re.match(r"^(chr\w+):(\d+)\s+([ACGTacgt]+)>([ACGTacgt]+)$", variant.strip())
        if not m:
            return {"enabled": True, "error": "expected chr21:25897620 C>T"}
        chrom, pos, ref, alt = m.group(1), int(m.group(2)), m.group(3).upper(), m.group(4).upper()
        adapter = AlphaGenomeAdapter()
        effects = adapter.predict(chrom, pos, ref, alt)
        effects.sort(key=lambda e: -abs(e.log2_fold_change))
        return {
            "enabled": True,
            "variant": f"{chrom}:{pos} {ref}>{alt}",
            "model": st["model"],
            "evidence": "predicted",
            "threshold_log2fc": 0.05,
            "scanned": adapter.last_scan,
            "effects": [
                {
                    "gene": e.gene,
                    "tissue": e.tissue,
                    "log2_fold_change": round(e.log2_fold_change, 4),
                    "direction": getattr(e.direction, "value", str(e.direction)),
                    "confidence": round(e.confidence, 3),
                }
                for e in effects
            ],
        }

    def lookup(self, variant: str) -> dict:
        """One variant through every layer: VEP anywhere, local trace when the chromosome is local."""
        from genomeos.genome import IndexedGenome, default_gencode
        from genomeos.genome.lookup import lookup, parse_variant

        try:
            chrom, pos, ref, alt = parse_variant(variant)
        except ValueError as e:
            raise ApiError(str(e)) from e
        ann = genome = None
        fa = self.root / "data" / "reference" / f"{chrom}.fa.gz"
        if default_gencode({chrom}) and fa.exists():
            ann = self._annotation_for(f"data/reference/{chrom}.fa.gz", chrom)
            genome = IndexedGenome(fa)
        try:
            return lookup(chrom, pos, ref, alt, ann, genome)
        finally:
            if genome:
                genome.close()

    def graph(self, gene: str, max_nodes: int = 40) -> dict:
        """One protein's neighbourhood in the local knowledge graph (cached definitions only)."""
        from genomeos.molecules.graph import cached_build

        if not gene or not gene.replace("-", "").replace("_", "").isalnum():
            raise ApiError("gene symbol required")
        g = cached_build()
        n = g.neighbourhood(gene.upper(), max_nodes)
        rels = ("associates", "member_of", "has_domain", "expressed_in")
        n["degree"] = {k: g.degree(gene.upper(), k) for k in rels}
        return n

    def rna(self, gene: str, chrom: str | None) -> dict:
        """Transcripts (local models) and GTEx expression per tissue for one gene."""
        from genomeos.genome import IndexedGenome, default_gencode
        from genomeos.molecules.rna import rna_report

        if not gene or not gene.replace("-", "").replace("_", "").isalnum():
            raise ApiError("gene symbol required")
        ann = genome = None
        if chrom and default_gencode({chrom}):
            ann = self._annotation_for(f"data/reference/{chrom}.fa.gz", chrom)
            fa = self.root / "data" / "reference" / f"{chrom}.fa.gz"
            genome = IndexedGenome(fa) if fa.exists() else None
        try:
            return rna_report(gene.upper(), ann, genome)
        finally:
            if genome:
                genome.close()

    def protein_definition(self, gene: str, refresh: bool = False) -> dict:
        """The federated protein definition (cached locally after the first compile)."""
        from genomeos.molecules import compile_protein, states_from_definition

        if not gene or not gene.replace("-", "").replace("_", "").isalnum():
            raise ApiError("gene symbol required")
        d = compile_protein(gene.upper(), refresh=refresh)
        d = dict(d)
        d["states"] = [st.to_dict() for st in states_from_definition(d)][:40]
        return d

    def pathway_kinetic(
        self, query: str, model: str | None, knockout: str | None, hours: float = 100.0, top: int = 6
    ) -> dict:
        """A curated BioModels ODE model for a pathway (by name or id), run on the in-house engine, with an
        optional knockout compared on final levels and peaks; the top species' time series for a chart."""
        import re

        from genomeos.molecules import biomodels

        hits, used = [], query
        if not model:
            if not query:
                raise ApiError("give a pathway name to search BioModels with, or a BIOMD… id")
            hits, used = biomodels.search_pathway(query, limit=8)
            if not hits:
                raise ApiError(f"no curated BioModels entry matches {query!r}")
            model = hits[0]["id"]
        if not re.fullmatch(r"(BIOMD|MODEL)\d+", model):
            raise ApiError("model id must look like BIOMD0000000010")
        path = biomodels.fetch(model, self.root / "data" / "knowledge" / "biomodels")
        dt = 0.01 if hours <= 500 else 0.05
        base = biomodels.run(path, duration=hours, dt=dt)
        out: dict = {
            "hits": hits,
            "query_used": used,
            "model": {k: v for k, v in base.items() if k not in ("series", "times", "levels")},
            "levels": base["levels"],
            "times": base["times"],
        }
        shown = [s for s, _ in sorted(base["levels"].items(), key=lambda kv: -kv[1]["peak"])[:top]]
        out["series"] = {s: base["series"][s] for s in shown}
        if knockout:
            acc = None
            if not re.fullmatch(r"[A-NR-Z][0-9][A-Z0-9]{3}[0-9]([A-Z][A-Z0-9]{2}[0-9])?", knockout.upper()):
                from genomeos.molecules import compile_protein

                try:
                    d = compile_protein(knockout.upper(), sources={"uniprot"})
                    acc = ((d["sections"].get("identity") or {}).get("items") or {}).get("accession")
                except Exception:  # noqa: BLE001 - match by name only
                    acc = None
            else:
                acc = knockout.upper()
            ko = biomodels.run(path, duration=hours, dt=dt, knockout=[knockout], accessions={knockout: acc})
            changed = biomodels.compare(base, ko)
            out["knockout"] = {
                "term": knockout,
                "accession": acc,
                "held": ko["knocked_out"],
                "changed": changed,
                "series": {r["species"]: ko["series"][r["species"]] for r in changed[:top]},
            }
            for r in changed[:top]:
                out["series"].setdefault(r["species"], base["series"][r["species"]])
        return out

    def pathway(self, pathway_id: str, knockout: str | None) -> dict:
        """A Reactome pathway as a reachability graph, optionally with one protein removed."""
        import re

        from genomeos.molecules.reactome import PathwayModel, fetch_pathway

        if not re.fullmatch(r"R-HSA-\d+", pathway_id or ""):
            raise ApiError("pathway id must look like R-HSA-69541")
        model = PathwayModel.from_sbml(fetch_pathway(pathway_id))
        out = {
            "summary": model.summary(),
            "reactions": [
                {
                    "id": r.reactome_id or r.id,
                    "name": r.name,
                    "inputs": [model.species[x].name for x in r.inputs],
                    "outputs": [model.species[x].name for x in r.outputs],
                    "catalysts": [model.species[x].name for x in r.catalysts],
                    "inhibitors": [model.species[x].name for x in r.inhibitors],
                }
                for r in model.reactions
            ],
        }
        if knockout:
            acc = knockout.upper()
            if not re.fullmatch(r"[A-NR-Z][0-9][A-Z0-9]{3}[0-9]([A-Z][A-Z0-9]{2}[0-9])?", acc):
                from genomeos.molecules import compile_protein

                ident = compile_protein(acc, sources={"uniprot"})["sections"].get("identity", {}).get("items")
                if not ident:
                    raise ApiError(f"no reviewed UniProt entry for {knockout}")
                acc = ident["accession"]
            out["knockout"] = model.knockout(acc)
            out["knockout"]["symbol"] = knockout.upper()
        return out

    def flow(self, gene: str, chrom: str, variant: str | None) -> dict:
        """DNA → RNA → protein trace of a gene's canonical transcript, optionally with one base change."""
        from genomeos.flow import trace_gene
        from genomeos.genome import IndexedGenome, default_gencode

        if not gene or not gene.replace("-", "").replace("_", "").isalnum():
            raise ApiError("gene symbol required")
        if not chrom or not default_gencode({chrom}):
            raise ApiError(f"no local annotation for {chrom or '?'}; chr21 and chrM are local")
        ann = self._annotation_for(f"data/reference/{chrom}.fa.gz", chrom)
        fa = self.root / "data" / "reference" / f"{chrom}.fa.gz"
        if not fa.exists():
            raise ApiError(f"no local sequence for {chrom}")
        try:
            g = ann.gene(gene.upper())
        except KeyError as e:
            raise ApiError(f"{gene} not on {chrom}") from e
        module = ann.to_module("flow")
        genome = IndexedGenome(fa)
        try:
            tr = trace_gene(genome, g, module.entities[g.id].transcripts)
        finally:
            genome.close()
        if tr is None:
            raise ApiError(f"{gene} has no transcripts")
        out = tr.to_dict()
        out["gene_start"], out["gene_end"] = g.locus.start, g.locus.end
        out["transcripts"] = len(g.transcripts)
        from genomeos.genome.regulation import regulation_of
        from genomeos.genome.regulatory import load_ccres

        ccres = load_ccres(chrom)
        if ccres:
            fai = IndexedGenome(fa)
            try:
                length = fai.lengths[chrom]
            finally:
                fai.close()
            out["regulation"] = regulation_of(g.symbol, chrom, ccres, ann, length)
        # every local individual's variants inside the gene: the test human and imported genomes
        from genomeos.genome.individuals import rows_in, sources

        people = []
        for name, vcf, evidence in sources(chrom, self.root / "data" / "individuals"):
            inside, coding, positions = 0, [], []
            for f in rows_in(vcf, g.locus.start + 1, g.locus.end):
                inside += 1
                pos = int(f[1])
                gt = f[9].split(":")[0] if len(f) > 9 else ""
                if len(positions) < 2000:
                    positions.append([pos - 1, gt])
                if tr is not None and len(f[3]) == 1 and len(f[4]) == 1 and len(coding) < 30:
                    sub = tr.substitute(pos - 1, f[3], f[4])
                    if sub.get("region") == "CDS":
                        coding.append(
                            {
                                "pos": pos,
                                "ref": f[3],
                                "alt": f[4],
                                "genotype": gt,
                                "consequence": sub.get("consequence"),
                                "hgvs_p": sub.get("hgvs_p"),
                                "residue": sub.get("residue"),
                            }
                        )
            people.append(
                {
                    "name": name,
                    "variants_in_gene": inside,
                    "positions": positions,
                    "coding_snvs": coding,
                    "evidence": evidence,
                }
            )
        out["individuals"] = people
        if variant:
            try:
                pos, change = (
                    variant.replace(",", "").split(":")[-1].split(maxsplit=1)
                    if " " in variant
                    else (variant.split(":")[-1][:-3], variant[-3:])
                )
                ref, alt = change.upper().split(">")
                out["variant"] = tr.substitute(int(pos) - 1, ref, alt)  # user gives 1-based
            except (ValueError, IndexError) as e:
                raise ApiError("variant format: POSITION REF>ALT, e.g. 25897620 G>A (1-based)") from e
        return out

    # ---- libraries -------------------------------------------------------

    def proteome_lib(self, gene: str | None = None) -> dict:
        """The packaged proteome library: its summary, or one protein's record and BioLang block."""
        from genomeos.lib.proteome import block, get, summary

        out: dict = {"summary": summary()}
        if gene:
            if not gene.replace("-", "").replace("_", "").isalnum():
                raise ApiError("gene symbol required")
            r = get(gene)
            out["protein"] = r
            out["block"] = block(gene) if r else None
        return out

    def libs(self) -> dict:
        return {
            "libraries": [
                {
                    "id": l.id,
                    "layer": l.layer,
                    "purpose": l.purpose,
                    "genes": list(l.genes),
                    "scale": l.scale,
                    "source": l.source,
                    "note": l.note,
                }
                for l in LIBRARIES.values()
            ]
        }


class Handler(BaseHTTPRequestHandler):
    api: Api  # set by serve()

    def log_message(self, fmt: str, *args) -> None:  # quieter console
        if self.server.verbose:  # type: ignore[attr-defined]
            super().log_message(fmt, *args)

    def _json(self, payload: dict, status: int = 200) -> None:
        body = json.dumps(payload).encode()
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _static(self, name: str) -> None:
        p = (STATIC / name).resolve()
        if STATIC.resolve() not in p.parents or not p.is_file():
            self._json({"error": "not found"}, 404)
            return
        body = p.read_bytes()
        ctype = {
            ".html": "text/html; charset=utf-8",
            ".js": "application/javascript; charset=utf-8",
            ".css": "text/css; charset=utf-8",
        }.get(p.suffix, "application/octet-stream")
        self.send_response(200)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _q(self, qs: dict, key: str, default=None):
        return qs.get(key, [default])[0]

    def do_GET(self) -> None:  # noqa: N802
        u = urlparse(self.path)
        qs = parse_qs(u.query)
        try:
            if u.path in ("/", "/index.html"):
                return self._static("index.html")
            if u.path.startswith("/static/") and "/" not in u.path[8:]:
                return self._static(u.path[8:])
            if u.path == "/api/files":
                return self._json(self.api.files())
            if u.path == "/api/module":
                return self._json(self.api.module_source(self._q(qs, "path", "")))
            if u.path == "/api/genome/info":
                return self._json(self.api.genome_info(self._q(qs, "path", "")))
            if u.path == "/api/genome/fetch":
                return self._json(self.api.genome_fetch(self._q(qs, "path", ""), self._q(qs, "locus", "")))
            if u.path == "/api/genome/orfs":
                return self._json(
                    self.api.genome_orfs(
                        self._q(qs, "path", ""),
                        int(self._q(qs, "min_aa", 50)),
                        self._q(qs, "table", "standard"),
                    )
                )
            if u.path == "/api/libs":
                return self._json(self.api.libs())
            if u.path == "/api/proteome_lib":
                return self._json(self.api.proteome_lib(self._q(qs, "gene")))
            if u.path == "/api/twins":
                return self._json(self.api.twins())
            if u.path == "/api/individuals":
                return self._json(self.api.individuals())
            if u.path == "/api/individual/screen":
                return self._json(self.api.individual_screen(self._q(qs, "name") or "HG002"))
            if u.path == "/api/individual/check":
                return self._json(
                    self.api.individual_check(self._q(qs, "name") or "HG002", self._q(qs, "chrom") or "chr21")
                )
            if u.path == "/api/individual/predict":
                return self._json(
                    self.api.individual_predict(
                        self._q(qs, "name") or "HG002",
                        self._q(qs, "gene") or "",
                        self._q(qs, "chrom") or "chr21",
                        int(self._q(qs, "max") or 30),
                    )
                )
            if u.path == "/api/individual/genes":
                return self._json(
                    self.api.individual_genes(
                        self._q(qs, "name") or "HG002",
                        self._q(qs, "chrom") or "chr21",
                        int(self._q(qs, "top") or 40),
                    )
                )
            if u.path == "/api/results":
                return self._json(self.api.results())
            if u.path == "/api/blocks":
                return self._json(
                    self.api.blocks(
                        self._q(qs, "path", ""),
                        self._q(qs, "chrom", ""),
                        int(self._q(qs, "start", 0)),
                        int(self._q(qs, "end", 0)),
                    )
                )
            if u.path == "/api/report":
                return self._json(self.api.report(self._q(qs, "gene", ""), self._q(qs, "chrom", "chr21")))
            if u.path == "/api/lookup":
                return self._json(self.api.lookup(self._q(qs, "variant", "")))
            if u.path == "/api/graph":
                return self._json(self.api.graph(self._q(qs, "gene", ""), int(self._q(qs, "max", 40))))
            if u.path == "/api/rna":
                return self._json(self.api.rna(self._q(qs, "gene", ""), self._q(qs, "chrom")))
            if u.path == "/api/protein_definition":
                return self._json(
                    self.api.protein_definition(self._q(qs, "gene", ""), self._q(qs, "refresh", "") == "1")
                )
            if u.path == "/api/pathway":
                return self._json(self.api.pathway(self._q(qs, "id", ""), self._q(qs, "knockout")))
            if u.path == "/api/pathway_kinetic":
                return self._json(
                    self.api.pathway_kinetic(
                        self._q(qs, "query", ""),
                        self._q(qs, "model"),
                        self._q(qs, "knockout"),
                        float(self._q(qs, "hours") or 100),
                    )
                )
            if u.path == "/api/flow":
                return self._json(
                    self.api.flow(
                        self._q(qs, "gene", ""), self._q(qs, "chrom", "chr21"), self._q(qs, "variant")
                    )
                )
            if u.path == "/api/protein":
                return self._json(self.api.protein(self._q(qs, "gene", ""), self._q(qs, "chrom")))
            if u.path == "/api/cancer":
                return self._json(self.api.cancer_knowledge(self._q(qs, "gene"), int(self._q(qs, "top", 30))))
            if u.path == "/api/jobs":
                return self._json(self.api.jobs())
            if u.path == "/api/progress":
                return self._json(self.api.progress())
            if u.path == "/api/features":
                return self._json(self.api.features())
            if u.path == "/api/predict":
                return self._json(self.api.predict(self._q(qs, "variant", "")))
            if u.path == "/api/anatomy/genome":
                return self._json(self.api.genome_wide())
            if u.path == "/api/unknown/genome":
                return self._json(self.api.unknown_wide())
            if u.path == "/api/individual":
                return self._json(self.api.individual_wide())
            if u.path == "/api/proteome":
                return self._json(self.api.proteome_summary())
            if u.path == "/api/anatomy":
                return self._json(self.api.anatomy(self._q(qs, "path"), self._q(qs, "chrom")))
            if u.path == "/api/grow":
                return self._json(
                    self.api.grow(
                        self._q(qs, "module", "data/organisms/celegans/embryo.bio"),
                        self._q(qs, "until", "800"),
                        int(self._q(qs, "depth", 3)),
                    )
                )
            if u.path == "/api/develop":
                return self._json(
                    self.api.develop(
                        self._q(qs, "model", "flag"),
                        int(self._q(qs, "width", 60)),
                        float(self._q(qs, "hours", 100)),
                    )
                )
            if u.path == "/api/age":
                kw = {k: v[0] for k, v in qs.items()}
                return self._json(self.api.age(**kw))
            return self._json({"error": "not found"}, 404)
        except ApiError as e:
            return self._json({"error": str(e)}, e.status)
        except Exception as e:  # noqa: BLE001
            return self._json({"error": f"{type(e).__name__}: {e}"}, 500)

    def do_POST(self) -> None:  # noqa: N802
        u = urlparse(self.path)
        length = int(self.headers.get("Content-Length", "0"))
        try:
            body = json.loads(self.rfile.read(length) or b"{}")
            if u.path == "/api/compile":
                return self._json(self.api.compile(body.get("source", "")))
            if u.path == "/api/cancer/tumour":
                return self._json(self.api.cancer_tumour(body.get("vcf", ""), int(body.get("deep", 5))))
            if u.path == "/api/therapeutics":
                return self._json(
                    self.api.therapeutics(
                        body.get("vcf", ""),
                        body.get("hla", ""),
                        body.get("rna", ""),
                        body.get("normal", ""),
                        int(body.get("top", 8)),
                        bool(body.get("offline", False)),
                        body.get("cohort", ""),
                    )
                )
            if u.path == "/api/cancer/compare":
                return self._json(
                    self.api.cancer_compare(
                        body.get("normal", ""),
                        body.get("tumour", ""),
                        body.get("genome", "data/reference/chr21.fa.gz"),
                        body.get("chrom", "chr21"),
                    )
                )
            if u.path == "/api/blocks/check":
                return self._json(
                    self.api.blocks_check_move(
                        body.get("blocks", []),
                        body.get("id", ""),
                        body.get("new_start", 0),
                        body.get("chrom_length"),
                    )
                )
            if u.path == "/api/jobs/heal":
                from genomeos import jobs

                return self._json(
                    jobs.heal(body.get("name", ""), self.api.root, force=bool(body.get("force"))).to_dict()
                )
            if u.path == "/api/jobs/start":
                return self._json(self.api.job_start(body.get("name", "")))
            if u.path == "/api/debug":
                return self._json(
                    self.api.debug(
                        body.get("source", ""),
                        body.get("breakpoints", []),
                        body.get("initial"),
                        body.get("until", 100),
                        body.get("explain"),
                        body.get("cell_type"),
                    )
                )
            if u.path == "/api/twin/measure":
                return self._json(self.api.twin_measure(body.get("name", "")))
            if u.path == "/api/individual/import":
                return self._json(
                    self.api.individual_import(
                        body.get("path", ""),
                        body.get("name", ""),
                        body.get("note", ""),
                        bool(body.get("replace")),
                    )
                )
            if u.path == "/api/twin/new":
                return self._json(
                    self.api.twin_new(
                        body.get("name", ""),
                        body.get("age", 0),
                        body.get("telomere"),
                        body.get("epigenetic_age"),
                        body.get("sex", "unknown"),
                    )
                )
            if u.path == "/api/twin/fork":
                return self._json(
                    self.api.twin_fork(
                        body.get("name", ""),
                        body.get("new_name", ""),
                        body.get("variants", []),
                        body.get("mutagen"),
                        body.get("pace"),
                        body.get("note", ""),
                    )
                )
            if u.path == "/api/twin/run":
                return self._json(
                    self.api.twin_run(
                        body.get("names", []),
                        body.get("cell_type", "fibroblast"),
                        body.get("years", 40),
                        body.get("cells", 300),
                        body.get("seed", 1),
                    )
                )
            if u.path == "/api/run":
                return self._json(
                    self.api.run(
                        body.get("source", ""),
                        body.get("hours", 50),
                        body.get("dt", 0.01),
                        body.get("initial"),
                        body.get("context"),
                        body.get("seed"),
                    )
                )
            return self._json({"error": "not found"}, 404)
        except ApiError as e:
            return self._json({"error": str(e)}, e.status)
        except Exception as e:  # noqa: BLE001
            return self._json({"error": f"{type(e).__name__}: {e}"}, 500)


def make_server(
    root: Path, host: str = "127.0.0.1", port: int = 8765, verbose: bool = False
) -> ThreadingHTTPServer:
    Handler.api = Api(root)
    # the supervisor: every minute, restart an auto-heal job that is not complete and is dead or stalled
    import threading

    from genomeos import jobs

    threading.Thread(target=jobs.supervise, args=(root, 60), daemon=True, name="job-supervisor").start()
    srv = ThreadingHTTPServer((host, port), Handler)
    srv.verbose = verbose  # type: ignore[attr-defined]
    return srv


def serve(
    root: Path | None = None,
    host: str = "127.0.0.1",
    port: int = 8765,
    open_browser: bool = False,
    verbose: bool = False,
) -> None:
    srv = make_server(root or Path.cwd(), host, port, verbose)
    url = f"http://{host}:{port}/"
    print(f"GenomeOS {__version__} web UI at {url}  (Ctrl-C to stop)")
    if open_browser:
        threading.Timer(0.5, lambda: webbrowser.open(url)).start()
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        srv.server_close()
