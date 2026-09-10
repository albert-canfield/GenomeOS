"""GenomeOS command line.

genomeos info    FASTA                 chromosome summary
genomeos orfs    FASTA [--min-aa N]    open reading frames -> proteins
genomeos compile MODULE.bio [-o OUT]   BioLang -> BioIR JSON
genomeos check   MODULE.bio            confidence / unknowns report
genomeos run     MODULE.bio [--hours]  execute a regulatory network
genomeos age     [--cell-type T]       cell ageing simulation
genomeos libs    [--layer L] [-v]      the biological libraries in the genome
genomeos serve   [--port 8765] [--open] light web UI on localhost
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from genomeos import __version__
from genomeos.genome import Genome, iter_fasta
from genomeos.ir import Module
from genomeos.lang import parse_file
from genomeos.lib import LAYERS, LIBRARIES, by_layer
from genomeos.runtime import (
    CELL_TYPES,
    STANDARD_CODE,
    VERTEBRATE_MITOCHONDRIAL_CODE,
    CellRuntime,
    Environment,
    NetworkRuntime,
    find_orfs,
)

CODON_TABLES = {"standard": STANDARD_CODE, "mito": VERTEBRATE_MITOCHONDRIAL_CODE}


def _table(rows: list[dict], cols: list[str]) -> str:
    widths = {c: max(len(c), *(len(str(r.get(c, ""))) for r in rows)) for c in cols}
    line = "  ".join(c.ljust(widths[c]) for c in cols)
    out = [line, "  ".join("-" * widths[c] for c in cols)]
    for r in rows:
        out.append("  ".join(str(r.get(c, "")).ljust(widths[c]) for c in cols))
    return "\n".join(out)


def _spark(values: list[float], width: int = 60) -> str:
    bars = "▁▂▃▄▅▆▇█"
    if not values:
        return ""
    step = max(1, len(values) // width)
    sample = values[::step][:width]
    lo, hi = min(sample), max(sample)
    if hi - lo < 1e-12:
        return bars[0] * len(sample)
    return "".join(bars[int((v - lo) / (hi - lo) * (len(bars) - 1))] for v in sample)


def cmd_info(args: argparse.Namespace) -> int:
    rows = []
    total = 0
    for name, _desc, seq in iter_fasta(args.fasta):
        fwd, rev = seq.telomeric_repeats()
        rows.append(
            {
                "chromosome": name,
                "length_bp": f"{len(seq):,}",
                "gc": f"{seq.gc_content():.3f}",
                "N": f"{seq.n_fraction():.3f}",
                "TTAGGG(+/-)": f"{fwd}/{rev}",
            }
        )
        total += len(seq)
    print(_table(rows, ["chromosome", "length_bp", "gc", "N", "TTAGGG(+/-)"]))
    print(f"\ntotal: {total:,} bp in {len(rows)} sequence(s)")
    return 0


def cmd_orfs(args: argparse.Namespace) -> int:
    table = CODON_TABLES[args.table]
    for name, _desc, seq in iter_fasta(args.fasta):
        for orf in find_orfs(seq, chrom=name, min_aa=args.min_aa, table=table):
            print(f"{orf.locus}\t{orf.length_aa} aa\t{orf.protein[:60]}{'...' if orf.length_aa > 60 else ''}")
    return 0


def _load_module(path: str) -> Module:
    p = Path(path)
    if p.suffix == ".json":
        return Module.from_dict(json.loads(p.read_text()))
    return parse_file(p)


def cmd_compile(args: argparse.Namespace) -> int:
    module = _load_module(args.module)
    text = json.dumps(module.to_dict(), indent=2)
    if args.output:
        Path(args.output).write_text(text)
        print(f"wrote {args.output}: {len(module.entities)} entities, {len(module.rules)} rules")
    else:
        print(text)
    return 0


def cmd_check(args: argparse.Namespace) -> int:
    module = _load_module(args.module)
    print(f"module {module.name}")
    print(
        f"  entities: {len(module.entities)}   rules: {len(module.rules)}   "
        f"parameters: {len(module.parameters)}   events: {len(module.events)}   "
        f"cell types: {len(module.cell_types())}"
    )
    if args.context:
        context = dict(kv.split("=", 1) for kv in args.context)
        active = module.active_rules(context)
        silenced = sorted(module.silenced_genes(context))
        print(
            f"  context {context}: {len(active)}/{len(module.rules)} rules active; "
            f"silenced: {silenced or 'none'}"
        )
        for r in active:
            print(f"    {r.id}")
    rep = module.confidence_report()
    print("  mean confidence by kind:")
    for kind, val in sorted(rep.items()):
        bar = "█" * int(val * 20)
        print(f"    {kind:<10} {bar:<20} {val:.2f}")
    unknowns = module.unknowns()
    if unknowns:
        print(f"  UNKNOWN regions: {len(unknowns)}")
        for u in unknowns:
            print(f"    {u.id}  {u.locus}")
    weak = [r for r in module.rules if r.confidence < 0.5]
    if weak:
        print(f"  rules with confidence < 0.5: {len(weak)}")
        for r in weak:
            print(f"    {r.id}  ({r.evidence.kind.value}: {r.evidence.source or 'no source'})")
    return 0


def _run_sbml(args: argparse.Namespace) -> int:
    from genomeos.runtime.sbml import SbmlModel, SbmlRuntime

    model = SbmlModel.from_file(args.module)
    rt = SbmlRuntime(model)
    traj = rt.run(duration=args.hours, dt=args.dt, record_every=max(1, int(args.hours / args.dt / 600)))
    print(
        f"SBML {model.id} {model.name!r}: {len(model.species)} species, {len(model.reactions)} reactions, "
        f"{len(model.assignments)} assignment rules; {args.hours} time units simulated"
    )
    for sp in traj.species:
        xs = traj.levels[sp]
        print(f"  {sp:<16} {_spark(xs)}  final={xs[-1]:10.3f}  peaks={traj.peaks(sp)}")
    return 0


def _run_boolean(args: argparse.Namespace) -> int:
    from genomeos.runtime.boolean import BooleanNetwork

    net = BooleanNetwork.from_file(args.module)
    init = {}
    for kv in args.init or []:
        k, v = kv.split("=", 1)
        init[k] = v.strip().lower() in ("1", "true", "on")
    print(f"Boolean network {args.module}: {len(net.nodes)} nodes")
    if init:
        att = net.attractor(init)
        kind = "fixed point" if len(att) == 1 else f"cycle of length {len(att)}"
        print(f"  from {init}: {kind}")
        for st in att:
            print("    " + " ".join(f"{n}={int(v)}" for n, v in st.items()))
    else:
        atts = net.attractors(samples=300, seed=args.seed or 0)
        print(f"  {len(atts)} distinct synchronous attractors from 300 random starts")
        for att in atts:
            on = sorted(n for n, v in att[0].items() if v)
            print(f"    length {len(att)}: first state on={on}")
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    if args.module.endswith(".xml"):
        return _run_sbml(args)
    if args.module.endswith(".bnet"):
        return _run_boolean(args)
    module = _load_module(args.module)
    context = dict(kv.split("=", 1) for kv in args.context) if args.context else {}
    vm = NetworkRuntime(module, context=context, seed=args.seed)
    initial = {}
    for kv in args.init or []:
        k, v = kv.split("=", 1)
        initial[k] = float(v)
    traj = vm.run(hours=args.hours, dt=args.dt, initial=initial)
    print(
        f"module {module.name}: {len(vm.active_rules)}/{len(module.rules)} rules active "
        f"in context {context or '{}'}; {args.hours} h simulated"
    )
    for s in traj.species:
        xs = traj.levels[s]
        print(f"  {s:<16} {_spark(xs)}  final={xs[-1]:8.2f}  peaks={traj.peaks(s)}")
    from genomeos.runtime.uncertainty import report_for_network

    print(report_for_network(module, vm.active_rules).format())
    if args.csv:
        with open(args.csv, "w") as fh:
            fh.write("time," + ",".join(traj.species) + "\n")
            for i, t in enumerate(traj.times):
                fh.write(f"{t:.4f}," + ",".join(f"{traj.levels[s][i]:.6g}" for s in traj.species) + "\n")
        print(f"  wrote {args.csv}")
    return 0


def cmd_age(args: argparse.Namespace) -> int:
    env = Environment(
        proliferation_factor=args.proliferation, mutagen_factor=args.mutagen, epigenetic_pace=args.pace
    )
    rt = CellRuntime(seed=args.seed, env=env)
    reports = rt.simulate_tissue(args.cell_type, args.years, n_cells=args.cells)
    rows = [
        {
            "age": f"{r.age_years:5.1f}",
            "senescent": f"{r.senescent_fraction:6.1%}",
            "telomere_bp": f"{r.mean_telomere_bp:8.0f}",
            "mutations": f"{r.mean_mutations:7.0f}",
            "epigenetic_age": f"{r.mean_epigenetic_age:6.1f}",
        }
        for r in reports
    ]
    print(f"{args.cell_type}: {args.cells} cells, {args.years} years, seed={args.seed}")
    print(_table(rows, ["age", "senescent", "telomere_bp", "mutations", "epigenetic_age"]))
    from genomeos.runtime.uncertainty import report_for_ageing

    print(report_for_ageing(rt.evidence_table(args.cell_type)).format())
    if args.evidence:
        print("\nparameters and evidence:")
        for p in rt.evidence_table(args.cell_type):
            print(
                f"  {p.name:<32} {p.value:>9g} {p.unit:<6} conf={p.confidence:.1f}  "
                f"[{p.evidence.kind.value}] {p.evidence.source}"
            )
    return 0


def cmd_annotate(args: argparse.Namespace) -> int:
    from genomeos.genome import Annotation

    chroms = set(args.chrom) if args.chrom else None
    ann = Annotation.from_gff3(args.gff3, chroms)
    print(f"{args.gff3}: {len(ann.genes):,} genes on {chroms or 'all chromosomes'}")
    rows = [{"gene_type": k, "count": v} for k, v in ann.summary().items()]
    print(_table(rows, ["gene_type", "count"]))
    if args.output:
        module = ann.to_module(args.name or f"gencode.{'_'.join(sorted(chroms)) if chroms else 'all'}")
        Path(args.output).write_text(json.dumps(module.to_dict()))
        print(f"wrote {args.output}: {len(module.entities):,} entities (BioIR v0.1)")
    return 0


def cmd_gene(args: argparse.Namespace) -> int:
    from genomeos.genome import Annotation
    from genomeos.runtime import coding_sequence, translate

    ann = Annotation.from_gff3(args.gff3, set(args.chrom) if args.chrom else None)
    try:
        g = ann.gene(args.symbol)
    except KeyError:
        print(f"gene {args.symbol!r} not found (restrict --chrom to speed up, or check the symbol)")
        return 1
    print(f"{g.symbol}  {g.id}  {g.type}  {g.locus}  {g.locus.length:,} bp  {len(g.transcripts)} transcripts")
    genome = Genome.from_fasta(args.genome) if args.genome else None
    table = CODON_TABLES[args.table]
    m = ann.to_module("gene")
    for tx in m.entities[g.id].transcripts:
        ttype = tx.attrs["transcript_type"]
        line = f"  {tx.attrs['name']:<16} {ttype:<24} exons={len(tx.exons):<3} cds={tx.cds or '-'}"
        if genome and tx.cds_segments:
            prot = translate(coding_sequence(genome, tx), table=table, initiator=True)
            line += f"  {len(prot)} aa  {prot[:40]}{'…' if len(prot) > 40 else ''}"
        print(line)
    return 0


def cmd_index(args: argparse.Namespace) -> int:
    from genomeos.genome import write_fai

    fai = write_fai(args.fasta)
    print(f"wrote {fai}")
    for line in Path(fai).read_text().splitlines()[:30]:
        name, length = line.split("\t")[:2]
        print(f"  {name:<12} {int(length):>14,} bp")
    return 0


def cmd_twin_build(args: argparse.Namespace) -> int:
    from genomeos.genome import IndexedGenome, Locus, iter_vcf, write_haplotypes

    genome = IndexedGenome(args.genome)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for chrom in args.chrom:
        variants = list(iter_vcf(args.vcf, {chrom}))
        ref = genome.fetch(Locus(chrom, 0, genome.lengths[chrom]))
        out = out_dir / f"{args.sample}_{chrom}.fa"
        stats = write_haplotypes(ref, chrom, variants, out, args.sample)
        phased = sum(1 for v in variants if v.phased)
        print(f"{chrom}: {len(variants):,} PASS variants ({phased:,} phased) -> {out}")
        for hap, st in stats.items():
            print(
                f"  {hap}: applied={st['applied']:,} snv={st['snv']:,} indel={st['indel']:,} "
                f"length_delta={st['length_delta']:+,} skipped_overlap={st['skipped_overlap']} "
                f"ref_mismatch={st['skipped_ref_mismatch']}"
            )
    return 0


def cmd_variant(args: argparse.Namespace) -> int:
    from genomeos.genome import Annotation, IndexedGenome, Variant
    from genomeos.runtime import classify_all

    ann = Annotation.from_gff3(args.gff3, {args.chrom})
    m = ann.to_module("v")
    genome = IndexedGenome(args.genome)
    pos = args.pos - 1 if args.one_based else args.pos
    v = Variant(args.chrom, pos, args.ref, (args.alt,), gt=(1, 1))
    hits = []
    for g in ann.protein_coding():
        if not (g.locus.start - 5000 <= pos < g.locus.end + 5000):
            continue
        txs = [t for t in m.entities[g.id].transcripts if t.cds_segments]
        for e in classify_all(genome, txs, v, CODON_TABLES[args.table]):
            hits.append((g.symbol, e))
    genome.close()
    if not hits:
        print(f"{args.chrom}:{pos} {args.ref}>{args.alt}: no protein-coding gene within 5 kb")
        return 0
    print(f"{args.chrom}:{pos} {args.ref}>{args.alt}  (0-based; most severe first)")
    seen = set()
    for sym, e in sorted(
        hits, key=lambda x: __import__("genomeos.runtime", fromlist=["severity"]).severity(x[1].consequence)
    ):
        key = (sym, e.transcript_id)
        if key in seen:
            continue
        seen.add(key)
        change = f"  {e.protein_change}" if e.protein_change else ""
        print(f"  {sym:<10} {e.transcript_id:<20} {e.consequence:<28}{change}  {e.note}")
    return 0


def cmd_cells(args: argparse.Namespace) -> int:
    from genomeos.knowledge import CellTypes

    cells = CellTypes.from_obo(args.obo)
    if args.search:
        for t in cells.search(args.search, args.limit):
            print(f"  {t.id}  {t.name}")
        return 0
    if args.lineage:
        term = cells.search(args.lineage, 1)
        if not term:
            print("no such cell type")
            return 1
        for line in cells.lineage(term[0].id):
            print("  " + line)
        return 0
    print(f"{args.obo}: {len(cells):,} cell types")
    return 0


def cmd_lr(args: argparse.Namespace) -> int:
    from genomeos.knowledge import LigandReceptorTable

    t = LigandReceptorTable.from_cellphonedb(args.interactions, args.genes)
    if args.gene:
        partners = t.partners(args.gene)
        print(f"{args.gene}: {len(partners)} partners")
        for p, cls_ in sorted(partners):
            print(f"  {p:<12} {cls_}")
        return 0
    m = t.to_module()
    print(f"{len(t):,} ligand-receptor pairs -> {len(m.rules):,} binds rules, {len(m.entities):,} proteins")
    if args.output:
        Path(args.output).write_text(json.dumps(m.to_dict()))
        print(f"wrote {args.output}")
    return 0


def _twin_path(name: str) -> Path:
    p = Path(name)
    return p if p.suffix == ".json" else Path("data/twins") / f"{name}.json"


def cmd_twin_new(args: argparse.Namespace) -> int:
    from genomeos.twin import Twin
    from genomeos.twin.twin import MeasuredState

    t = Twin(
        args.name,
        genome=args.genome or "",
        vcf=args.vcf or "",
        sex=args.sex,
        measured=MeasuredState(args.age, args.telomere, args.epigenetic_age),
    )
    p = t.save()
    print(f"created twin {t.name} -> {p}")
    return 0


def cmd_twin_fork(args: argparse.Namespace) -> int:
    from genomeos.twin import Twin

    t = Twin.load(_twin_path(args.name)).fork(args.new_name, args.note or "")
    if args.variant:
        for spec in args.variant:
            gene, _, cons = spec.partition(":")
            applied = t.add_variant(gene, cons or "nonsense")
            mods = applied or "none (not a timer/maintenance LoF)"
            print(f"  variant {gene} {cons or 'nonsense'} -> modifiers {mods}")
    if args.mutagen is not None:
        t.environment.mutagen_factor = args.mutagen
    if args.pace is not None:
        t.environment.epigenetic_pace = args.pace
    p = t.save()
    print(f"forked {args.name} -> {t.name} ({p})")
    return 0


def cmd_twin_run(args: argparse.Namespace) -> int:
    from genomeos.twin import Twin, diff_runs

    runs = []
    for name in args.names:
        t = Twin.load(_twin_path(name))
        run = t.run(args.cell_type, years=args.years, cells=args.cells, seed=args.seed)
        runs.append(run)
        print(
            f"twin {t.name}  ({args.cell_type}, {args.years} y from age {t.measured.chronological_age})"
            f"{'  modifiers=' + str(run.modifiers) if run.modifiers else ''}"
        )
        rows = [
            {
                "age": f"{r.age_years:5.1f}",
                "senescent": f"{r.senescent_fraction:6.1%}",
                "telomere_bp": f"{r.mean_telomere_bp:8.0f}",
                "mutations": f"{r.mean_mutations:7.0f}",
                "epigenetic_age": f"{r.mean_epigenetic_age:6.1f}",
            }
            for r in run.reports[:: max(1, len(run.reports) // 8)]
        ]
        print(_table(rows, ["age", "senescent", "telomere_bp", "mutations", "epigenetic_age"]))
        print(run.uncertainty.format())
        if args.output:
            out = Path(args.output) / f"{t.name}.{args.cell_type}.json"
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text(json.dumps(run.to_dict()))
    if len(runs) == 2:
        print(f"\ndiff {runs[1].twin} - {runs[0].twin} at {runs[0].final.age_years:.0f} y:")
        for k, v in diff_runs(runs[0], runs[1]).items():
            print(f"  {k:<20} {v['a']:10.2f} -> {v['b']:10.2f}   delta {v['delta']:+.2f}")
    return 0


def cmd_clock(args: argparse.Namespace) -> int:
    from genomeos.twin.clocks import Clock, MethylationMatrix

    clocks = [Clock.horvath(), Clock.hannum()] if args.clock == "all" else [getattr(Clock, args.clock)()]
    path = args.methylation
    m = (
        MethylationMatrix.from_geo_series_matrix(path, cpgs={c for k in clocks for c in k.coefficients})
        if "series_matrix" in path
        else MethylationMatrix.from_csv(path)
    )
    samples = args.sample or m.samples[: args.limit]
    for s in samples:
        betas = m.sample(s)
        known = m.metadata.get(s, {}).get("age", "")
        parts = []
        for k in clocks:
            r = k.predict(betas)
            parts.append(
                f"{k.name}={r.age:5.1f}y (coverage {r.coverage:.0%}, conf {r.parameter.confidence:.2f})"
            )
        print(f"  {s:<12} {'age ' + known + '  ' if known else ''}{'  '.join(parts)}")
    return 0


def cmd_telomere(args: argparse.Namespace) -> int:
    from genomeos.genome.telomere import estimate_file

    est = estimate_file(args.reads, k=args.k, max_reads=args.max_reads)
    if args.save:
        from genomeos.results import save_result

        p = est.parameter
        path = save_result(
            args.save,
            {
                "source": args.reads,
                "reads": est.reads,
                "telomeric_reads": est.telomeric_reads,
                "mean_read_length": round(est.read_length, 1),
                "telomere_bp": round(est.telomere_bp),
                "fraction": est.fraction,
                "k": args.k,
                "confidence": p.confidence,
                "evidence": {"kind": p.evidence.kind.value, "source": p.evidence.source},
            },
        )
        print(f"  saved {path}")
    print(
        f"{args.reads}: {est.reads:,} reads, {est.telomeric_reads:,} telomeric ({est.fraction:.2e}), "
        f"mean read length {est.read_length:.0f}"
    )
    p = est.parameter
    print(f"  estimated mean telomere length: {p.value:,.0f} bp  conf={p.confidence:.1f}")
    print(f"  [{p.evidence.kind.value}] {p.evidence.source}")
    return 0


def cmd_debug(args: argparse.Namespace) -> int:
    from genomeos.runtime.debugger import AgeingDebugger, NetworkDebugger

    if args.target in CELL_TYPES:
        dbg = AgeingDebugger(args.target, seed=args.seed or 1)
        for b in args.breakpoint or []:
            print(f"breakpoint: {dbg.add_breakpoint(b)}")
        hit = dbg.step(years=args.until)
        print(f"stopped at {'breakpoint ' + str(hit) if hit else f'{args.until} years (no breakpoint hit)'}")
        print("state:", {k: round(v, 2) for k, v in dbg.session.state.items()})
        print("trace (last 8):")
        for line in dbg.session.trace[-8:]:
            print("  " + line.format())
        print("explain:")
        for line in dbg.explain():
            print("  " + line.format())
        return 0
    module = _load_module(args.target)
    context = dict(kv.split("=", 1) for kv in args.context) if args.context else {}
    dbg = NetworkDebugger(module, context=context, dt=args.dt)
    if args.init:
        dbg.set_initial({k: float(v) for k, v in (kv.split("=", 1) for kv in args.init)})
    for b in args.breakpoint or []:
        print(f"breakpoint: {dbg.add_breakpoint(b)}")
    hit = dbg.run_until_break(max_hours=args.until)
    print(
        f"stopped at t={dbg.session.time:.2f} h: {'breakpoint ' + str(hit) if hit else 'no breakpoint hit'}"
    )
    print("state:", {k: round(v, 3) for k, v in dbg.session.state.items()})
    for sp in args.explain or ([hit.variable] if hit else []):
        print(f"explain {sp}:")
        for line in dbg.explain(sp):
            print("  " + line.format())
    return 0


def cmd_develop(args: argparse.Namespace) -> int:
    if args.model == "flag":
        from genomeos.runtime.spatial import french_flag

        rt = french_flag(width=args.width, height=6, hours=args.hours)
        print("French flag (Wolpert 1969): morphogen source at the left edge")
        for row in rt.type_map()[:3]:
            print("  " + row)
        print("  bands:", rt.bands_along_x(), "census:", rt.census())
        for p in rt.parameters:
            print(f"  {p.name}={p.value} [{p.evidence.kind.value}] {p.evidence.source} conf={p.confidence}")
        return 0
    if args.model == "segmentation":
        from genomeos.runtime.segmentation import run_segmentation

        r = run_segmentation(length=args.width, hours=args.hours)
        print(f"segmentation clock: period {r.period_h} h, {r.frozen_cells} cells frozen")
        print(f"  {r.count} segments (expected {r.expected:.1f})")
        row = "".join(
            "|" if any(s == x for s, _ in r.segments) else "A" if ph < 3.14159 else "P"
            for x, ph in enumerate(r.frozen_phase)
            if ph == ph
        )
        print("  " + row)
        print(r.uncertainty().format())
        return 0
    from genomeos.runtime.gastrulation import run_gastrulation

    r = run_gastrulation(cells=args.width, hours=args.hours)
    print("gastrulation: NODAL source at the left, fates along the axis")
    print(
        "  "
        + "".join(f[0].upper() for f in r.fates)
        + "   (E=endoderm, M=mesoderm, E…=ectoderm shown as 'E'/'M'/'E')"
    )
    print("  proportions:", {k: round(v, 2) for k, v in r.proportions().items()})
    print(r.uncertainty().format())
    return 0


def cmd_forge(args: argparse.Namespace) -> int:
    from genomeos.forge import Knob, forge, period_of

    module = _load_module(args.module)
    initial = {k: float(v) for k, v in (kv.split("=", 1) for kv in (args.init or []))}
    knobs = []
    for spec in args.vary:
        path, _, rng = spec.partition("@")
        lo, _, hi = rng.partition("..")
        knobs.append(Knob(path, float(lo or 0.01), float(hi or 1000)))
    kind, _, target = args.target.partition("=")
    target_v = float(target)
    species = args.species

    def loss(traj):
        if kind == "period":
            p = period_of(traj, species)
            return abs(p - target_v) / target_v if p else 10.0
        if kind == "level":
            return abs(traj.final()[species] - target_v) / max(target_v, 1e-9)
        raise SystemExit("target must be period=<h> or level=<value>")

    def oscillates(traj):
        return traj.peaks(species) >= 3

    res = forge(
        module,
        knobs,
        loss,
        [oscillates] if kind == "period" else [],
        hours=args.hours,
        dt=args.dt,
        initial=initial,
        iterations=args.iterations,
        restarts=args.restarts,
        seed=args.seed,
    )
    status = "feasible" if res.best.feasible else "INFEASIBLE"
    print(f"BioForge: {res.evaluations} simulations; best loss {res.best.loss:.3f} ({status})")
    for path, before, after in res.changes():
        print(f"  {path:<36} {before:10.4g} -> {after:10.4g}")
    print("  every changed value now carries evidence 'predicted: BioForge search' at confidence <= 0.3")
    if args.output:
        Path(args.output).write_text(json.dumps(res.best.module.to_dict(), indent=1))
        print(f"  wrote {args.output}")
    return 0


def cmd_organism(args: argparse.Namespace) -> int:
    from genomeos.organism import run_lineage

    lin = run_lineage(until_min=args.minutes)
    alive = lin.count_at(args.minutes)
    print(f"C. elegans early lineage to {args.minutes:.0f} min: {len(lin.cells)} cells born, {alive} alive")
    for t in (0.5, 23, 50, 100, args.minutes):
        if t <= args.minutes:
            by = {lg: lin.lineage_count_at(lg, t) for lg in ("AB", "MS", "E", "C", "D", "P")}
            print(f"  t={t:5.0f} min  cells={lin.count_at(t):3d}  {by}")
    for line in lin.tree(depth=args.depth):
        print("  " + line)
    src = lin.module.parameters["ab_cycle_min"].evidence.source
    print(f"  cycle times: [experimental] {src}")
    return 0


def cmd_data(args: argparse.Namespace) -> int:
    from genomeos import storage
    from genomeos.results import list_results

    if args.data_cmd == "status":
        st = storage.status()
        print("data directories:")
        for name, size in st["dirs"].items():
            print(f"  {name:<12} {storage.human(size):>10}")
        print("largest files:")
        for path, size in st["largest"][:8]:
            print(f"  {storage.human(size):>10}  {path}")
        d = st["disk"]
        print(f"disk: {storage.human(d['free'])} free of {storage.human(d['total'])}")
        print(f"result summaries: {st['results']}")
        print("distillers:")
        for name, m in storage.manifest().items():
            print(f"  {name:<26} {'summary ok ' if m['summary_exists'] else 'no summary '} {m['describe']}")
        return 0
    if args.data_cmd == "distil":
        for name, what in storage.distil(args.only, force=args.force):
            print(f"  {name:<26} {what}")
        return 0
    if args.data_cmd == "clean":
        rows = storage.clean(dry_run=not args.yes)
        total = sum(size for _, size, _ in rows)
        for path, size, action in rows:
            print(f"  {action:<13} {storage.human(size):>10}  {path}")
        print(
            f"  {'would free' if not args.yes else 'freed'} {storage.human(total)}"
            + ("" if args.yes else "  (add --yes to delete)")
        )
        return 0
    for r in list_results():
        print(f"  {r['name']:<28} {r['date']}  {storage.human(r['size']):>8}  keys: {', '.join(r['keys'])}")
    return 0


def cmd_libs(args: argparse.Namespace) -> int:
    from genomeos.lib import KnowledgeBase

    kb = None
    if args.data or args.verify or args.gene:
        if not KnowledgeBase.available():
            print("knowledge files missing: run scripts/fetch_reference.sh knowledge")
            return 1
        kb = KnowledgeBase()
    if args.gene:
        libs = kb.libraries_of(args.gene)
        print(f"{args.gene}: member of {len(libs)} libraries by data")
        for lib_id in libs:
            print(f"  {lib_id}")
        terms = sorted(kb.annotations.direct.get(args.gene, ()))
        print(f"  direct GO annotations: {len(terms)}")
        for t in terms[:15]:
            print(f"    {t}  {kb.go.name(t)}")
        return 0
    if args.verify:
        bad = 0
        for v in kb.verify_all():
            flag = "ok " if v.ok else "MISS"
            bad += not v.ok
            extra = ""
            if v.missing:
                extra += f"  missing={list(v.missing)}"
            if v.unknown_terms:
                extra += f"  unknown_terms={list(v.unknown_terms)}"
            print(f"{flag} {v.library:<40} members={v.members:>5}  checked={len(v.checked):>2}{extra}")
        print(f"{bad} libraries with problems")
        return 1 if bad else 0
    layers = [args.layer] if args.layer else list(LAYERS)
    for layer in layers:
        libs = by_layer(layer)
        print(f"[{layer}]  {len(libs)} libraries")
        for lib in libs:
            count = (
                f"  [{len(kb.members(lib)):,} genes by data]" if kb and (lib.go_terms or lib.reactome) else ""
            )
            print(f"  {lib.id:<40} {lib.purpose}{count}")
            if args.verbose:
                if lib.genes:
                    print(f"      genes:  {', '.join(lib.genes)}")
                if lib.scale:
                    print(f"      scale:  {lib.scale}")
                if lib.source:
                    print(f"      source: {lib.source}")
                if lib.note:
                    print(f"      note:   {lib.note}")
        print()
    total_genes = len({g for l in LIBRARIES.values() for g in l.genes})
    print(f"{len(LIBRARIES)} libraries, {total_genes} representative genes")
    return 0


def cmd_serve(args: argparse.Namespace) -> int:
    from genomeos.web import serve

    serve(Path.cwd(), host=args.host, port=args.port, open_browser=args.open, verbose=args.verbose)
    return 0


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="genomeos", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--version", action="version", version=f"genomeos {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("info", help="summarise a FASTA file")
    p.add_argument("fasta")
    p.set_defaults(fn=cmd_info)

    p = sub.add_parser("orfs", help="find open reading frames")
    p.add_argument("fasta")
    p.add_argument("--min-aa", type=int, default=50)
    p.add_argument(
        "--table",
        choices=sorted(CODON_TABLES),
        default="standard",
        help="codon table: standard (nuclear) or mito (vertebrate mitochondrial)",
    )
    p.set_defaults(fn=cmd_orfs)

    p = sub.add_parser("compile", help="BioLang -> BioIR JSON")
    p.add_argument("module")
    p.add_argument("-o", "--output")
    p.set_defaults(fn=cmd_compile)

    p = sub.add_parser("check", help="confidence and unknowns report")
    p.add_argument("module")
    p.add_argument("--context", nargs="*", metavar="KEY=VALUE", help="e.g. cell_type=Neuron")
    p.set_defaults(fn=cmd_check)

    p = sub.add_parser("run", help="execute a regulatory network module")
    p.add_argument("module")
    p.add_argument("--hours", type=float, default=50.0)
    p.add_argument("--dt", type=float, default=0.01)
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--context", nargs="*", metavar="KEY=VALUE")
    p.add_argument("--init", nargs="*", metavar="SPECIES=LEVEL")
    p.add_argument("--csv")
    p.set_defaults(fn=cmd_run)

    p = sub.add_parser("annotate", help="compile a GFF3 annotation (GENCODE) to BioIR")
    p.add_argument("gff3")
    p.add_argument("--chrom", nargs="*", help="restrict to chromosomes, e.g. chr21 chrM")
    p.add_argument("--name", help="module name")
    p.add_argument("-o", "--output", help="write BioIR JSON")
    p.set_defaults(fn=cmd_annotate)

    p = sub.add_parser("gene", help="show a gene's transcripts and translated proteins")
    p.add_argument("symbol")
    p.add_argument("--gff3", required=True)
    p.add_argument("--genome", help="FASTA of the gene's chromosome")
    p.add_argument("--chrom", nargs="*")
    p.add_argument("--table", choices=sorted(CODON_TABLES), default="standard")
    p.set_defaults(fn=cmd_gene)

    p = sub.add_parser("index", help="build a .fai index for random access (decompresses .gz)")
    p.add_argument("fasta")
    p.set_defaults(fn=cmd_index)

    twin = sub.add_parser("twin", help="build and run an individual's genome").add_subparsers(
        dest="twin_cmd", required=True
    )
    p = twin.add_parser("new", help="create a twin with measured state")
    p.add_argument("name")
    p.add_argument("--genome")
    p.add_argument("--vcf")
    p.add_argument("--sex", default="unknown")
    p.add_argument("--age", type=float, default=0.0, help="chronological age in years")
    p.add_argument("--telomere", type=float, help="measured mean telomere length, bp")
    p.add_argument("--epigenetic-age", type=float, help="measured epigenetic age, years")
    p.set_defaults(fn=cmd_twin_new)

    p = twin.add_parser("fork", help="copy a twin with variants or environment changes")
    p.add_argument("name")
    p.add_argument("new_name")
    p.add_argument("--variant", nargs="*", metavar="GENE:CONSEQUENCE", help="e.g. TERT:nonsense")
    p.add_argument("--mutagen", type=float)
    p.add_argument("--pace", type=float)
    p.add_argument("--note")
    p.set_defaults(fn=cmd_twin_fork)

    p = twin.add_parser("run", help="run one or two twins forward and diff them")
    p.add_argument("names", nargs="+")
    p.add_argument("--cell-type", default="fibroblast", choices=sorted(CELL_TYPES))
    p.add_argument("--years", type=float, default=40.0)
    p.add_argument("--cells", type=int, default=500)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("-o", "--output", help="directory for run JSON")
    p.set_defaults(fn=cmd_twin_run)

    p = twin.add_parser("build", help="apply a VCF to a reference to produce hap1/hap2 FASTA")
    p.add_argument("--genome", required=True, help="reference FASTA (indexed on first use)")
    p.add_argument("--vcf", required=True)
    p.add_argument("--chrom", nargs="+", required=True)
    p.add_argument("--sample", default="HG002")
    p.add_argument("--out-dir", default="data/twins")
    p.set_defaults(fn=cmd_twin_build)

    p = sub.add_parser("variant", help="effect of a variant on the genes around it")
    p.add_argument("chrom")
    p.add_argument("pos", type=int)
    p.add_argument("ref")
    p.add_argument("alt")
    p.add_argument("--gff3", required=True)
    p.add_argument("--genome", required=True)
    p.add_argument("--one-based", action="store_true", help="POS is 1-based as in a VCF")
    p.add_argument("--table", choices=sorted(CODON_TABLES), default="standard")
    p.set_defaults(fn=cmd_variant)

    p = sub.add_parser("cells", help="cell types from the Cell Ontology")
    p.add_argument("--obo", default="data/knowledge/cl-basic.obo")
    p.add_argument("--search")
    p.add_argument("--lineage", help="print the is_a lineage of the best match")
    p.add_argument("--limit", type=int, default=20)
    p.set_defaults(fn=cmd_cells)

    p = sub.add_parser("lr", help="ligand-receptor pairs (CellPhoneDB) as binds rules")
    p.add_argument("--interactions", default="data/knowledge/cellphonedb_interaction_input.csv")
    p.add_argument("--genes", default="data/knowledge/cellphonedb_gene_input.csv")
    p.add_argument("--gene", help="list partners of a gene")
    p.add_argument("-o", "--output", help="write the BioIR module")
    p.set_defaults(fn=cmd_lr)

    p = sub.add_parser("clock", help="epigenetic age from a methylation matrix")
    p.add_argument("methylation", help="CSV (CpG rows × sample columns) or GEO series_matrix.txt(.gz)")
    p.add_argument("--clock", choices=["horvath", "hannum", "all"], default="all")
    p.add_argument("--sample", nargs="*")
    p.add_argument("--limit", type=int, default=10)
    p.set_defaults(fn=cmd_clock)

    p = sub.add_parser(
        "telomere", help="telomere length from reads: FASTQ/BAM on disk, or a FASTQ URL streamed"
    )
    p.add_argument("reads")
    p.add_argument("-k", type=int, default=7, help="TTAGGG repeats per read to call it telomeric")
    p.add_argument("--max-reads", type=int)
    p.add_argument("--save", metavar="NAME", help="write a result summary to data/results/NAME.json")
    p.set_defaults(fn=cmd_telomere)

    p = sub.add_parser("debug", help="step a module or a cell with breakpoints and an evidence trace")
    p.add_argument("target", help="a .bio module, or a cell type for the ageing debugger")
    p.add_argument(
        "-b", "--breakpoint", nargs="*", metavar="COND", help="e.g. 'TetR > 50' or 'divisions >= 10'"
    )
    p.add_argument("--until", type=float, default=100.0, help="max hours (module) or years (cell)")
    p.add_argument("--dt", type=float, default=0.01)
    p.add_argument("--init", nargs="*", metavar="SPECIES=LEVEL")
    p.add_argument("--context", nargs="*", metavar="KEY=VALUE")
    p.add_argument("--explain", nargs="*", metavar="SPECIES")
    p.add_argument("--seed", type=int)
    p.set_defaults(fn=cmd_debug)

    p = sub.add_parser("develop", help="developmental models: flag, segmentation, gastrulation")
    p.add_argument("model", choices=["flag", "segmentation", "gastrulation"])
    p.add_argument("--width", type=int, default=60)
    p.add_argument("--hours", type=float, default=100.0)
    p.set_defaults(fn=cmd_develop)

    p = sub.add_parser("forge", help="design under constraints: tune a module toward a target")
    p.add_argument("module")
    p.add_argument("--target", required=True, help="period=<hours> or level=<value>")
    p.add_argument("--species", required=True, help="species the target refers to, e.g. TetR")
    p.add_argument(
        "--vary",
        nargs="+",
        required=True,
        metavar="KNOB[@lo..hi]",
        help="gene:tetR.max, rule:<id>.strength, param:translation_rate",
    )
    p.add_argument("--init", nargs="*", metavar="SPECIES=LEVEL")
    p.add_argument("--hours", type=float, default=120.0)
    p.add_argument("--dt", type=float, default=0.02)
    p.add_argument("--iterations", type=int, default=40)
    p.add_argument("--restarts", type=int, default=2)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("-o", "--output", help="write the designed module as BioIR JSON")
    p.set_defaults(fn=cmd_forge)

    p = sub.add_parser("organism", help="minimal organism: C. elegans early lineage")
    p.add_argument("--minutes", type=float, default=150.0)
    p.add_argument("--depth", type=int, default=2)
    p.set_defaults(fn=cmd_organism)

    data = sub.add_parser("data", help="storage economics: status, distil, clean, results").add_subparsers(
        dest="data_cmd", required=True
    )
    data.add_parser("status", help="disk use and what can be distilled")
    p = data.add_parser("distil", help="turn raw downloads into result summaries")
    p.add_argument("--only", nargs="*")
    p.add_argument("--force", action="store_true", help="recompute existing summaries")
    p = data.add_parser("clean", help="delete raw inputs that already have a summary")
    p.add_argument("--yes", action="store_true", help="actually delete (default is a dry run)")
    data.add_parser("results", help="list result summaries")
    for sp in data.choices.values():
        sp.set_defaults(fn=cmd_data)

    p = sub.add_parser("libs", help="list the biological libraries found in the genome")
    p.add_argument("--layer", choices=LAYERS)
    p.add_argument("-v", "--verbose", action="store_true")
    p.add_argument("--data", action="store_true", help="count members from GO/Reactome data")
    p.add_argument("--verify", action="store_true", help="check catalogue genes against the data")
    p.add_argument("--gene", help="which libraries a gene belongs to, by data")
    p.set_defaults(fn=cmd_libs)

    p = sub.add_parser("serve", help="start the light web UI on localhost")
    p.add_argument("--host", default="127.0.0.1")
    p.add_argument("--port", type=int, default=8765)
    p.add_argument("--open", action="store_true", help="open the browser")
    p.add_argument("-v", "--verbose", action="store_true", help="log requests")
    p.set_defaults(fn=cmd_serve)

    p = sub.add_parser("age", help="simulate cell ageing clocks")
    p.add_argument("--cell-type", default="fibroblast", choices=sorted(CELL_TYPES))
    p.add_argument("--years", type=float, default=90.0)
    p.add_argument("--cells", type=int, default=500)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--proliferation", type=float, default=1.0)
    p.add_argument("--mutagen", type=float, default=1.0)
    p.add_argument("--pace", type=float, default=1.0)
    p.add_argument("--evidence", action="store_true", help="print parameter provenance")
    p.set_defaults(fn=cmd_age)
    return ap


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
