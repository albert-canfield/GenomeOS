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
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.

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
    widths = {c: max([len(c), *(len(str(r.get(c, ""))) for r in rows)]) for c in cols}  # rows may be empty
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

    if args.bai:
        from genomeos.genome.telomere import estimate_remote
        from genomeos.results import save_result

        chroms = args.chrom or [f"chr{i}" for i in range(1, 23)] + ["chrX"]
        r = estimate_remote(
            args.reads,
            args.bai,
            chroms,
            window=args.window,
            sample_bytes=int(args.sample_mb * 1e6),
            k=args.k,
            progress=lambda m: print("  " + m, flush=True),
        )
        p = r.pop("parameter")
        r["confidence"] = p.confidence
        r["evidence"] = {"kind": p.evidence.kind.value, "source": p.evidence.source}
        with_reads = [e for e in r["ends"] if e["reads"]]
        print(
            f"{args.reads}: {r['reads_total_from_index']:,} reads in the index "
            f"({r['reads_unmapped']:,} unmapped); "
            f"{r['unmapped_sampled']:,} unmapped sampled, {r['unmapped_sample_telomeric']:,} telomeric "
            f"({r['unmapped_telomeric_fraction']}); {len(with_reads)} chromosome ends read, "
            f"mean coverage {r['mean_end_coverage']}x; {r['bytes_fetched'] / 1e6:.0f} MB "
            f"in {r['requests']} requests"
        )
        print(f"  telomere ≈ {r['telomere_bp']:,} bp per end (TelSeq over ranges, confidence {p.confidence})")
        if args.save:
            print(f"  saved {save_result(args.save, r)}")
        return 0
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


def cmd_grow(args: argparse.Namespace) -> int:
    """Grow an organism from one cell: BioLang v0.3 organism program -> Body runtime -> report."""
    from genomeos.ir import to_minutes
    from genomeos.lang import parse_file
    from genomeos.organism import REFERENCES, ReferenceLineage
    from genomeos.organism.diff import compare
    from genomeos.runtime.body import Body

    module = parse_file(args.module)
    if module.organism is None:
        print(f"{args.module}: no organism block (see docs/BIOLANG-v0.3.md)", file=sys.stderr)
        return 2
    o = module.organism
    if args.design is not None:
        from genomeos.organism.forge import run_designs

        names = [x for x in args.design.split(",") if x] or None
        results = run_designs(module, names, seed=args.seed)
        if not results:
            print(f"{args.module}: no design blocks" + (f" named {names}" if names else ""))
            return 1
        for r in results:
            print(r.format())
        if args.json:
            Path(args.json).write_text(json.dumps([r.to_dict() for r in results], indent=1, default=str))
            print(f"  wrote {args.json}")
        return 0
    if args.experiments is not None:
        from genomeos.organism.experiment import run_experiments

        names = [x for x in args.experiments.split(",") if x] or None
        results = run_experiments(module, names, seed=args.seed)
        if not results:
            print(f"{args.module}: no experiment blocks" + (f" named {names}" if names else ""))
            return 1
        for r in results:
            print(r.format())
        if args.json:
            Path(args.json).write_text(json.dumps([r.to_dict() for r in results], indent=1, default=str))
            print(f"  wrote {args.json}")
        return 0
    parts = str(args.until).split()
    args.until = to_minutes(float(parts[0]), parts[1] if len(parts) > 1 else "min")
    body = Body(module, seed=args.seed, max_cells=args.max_cells, means=args.means).run(until=args.until)
    s = body.summary()
    when = f"{args.until:g} min"
    if args.until > 10_000:
        when = f"{args.until / 1440:.0f} d ({args.until / 525960:.1f} yr)"
    lost = ""
    if s["culled"]:
        lost = f", {s['culled']:.3g} lost to turnover, {s['turnover_per_day']:.3g} replaced per day"
    print(
        f"{o.name} ({o.species or 'species not stated'}): {s['cells_born']} cells born, "
        f"{s['alive']:.4g} alive at {when}, {s['deaths']} deaths{lost}; "
        f"{len(module.decisions)} decisions / {len(module.timers)} timers / {len(module.signals())} signals"
    )
    marks = (0.5, 50, 100, 200, 350, 500, 800, 2000, 4000, 5700)
    if args.until > 10_000:  # long runs: days, months, years
        days = (1, 5, 14, 21, 56, 266)
        years = (1, 5, 18, 40, 80)
        marks = tuple(to_minutes(x, "d") for x in days) + tuple(to_minutes(x, "yr") for x in years)
    checkpoints = [t for t in marks if t <= args.until]
    if args.until not in checkpoints:
        checkpoints.append(args.until)
    for t in checkpoints:
        parts = []
        for what in o.observe or ["count"]:
            if what == "count":
                parts.append(f"cells={body.count_at(t):g}")
            elif what == "deaths":
                parts.append(f"deaths={body.deaths_by(t)}")
            elif what == "fates":
                top = list(body.fates_at(t).items())[:4]
                parts.append("fates=" + ", ".join(f"{k} {v:g}" for k, v in top))
            elif what.startswith("lineage "):
                lg = what.split(None, 1)[1]
                parts.append(f"{lg}={body.lineage_count_at(lg, t):g}")
        label = f"t={t:7.1f} min" if args.until <= 10_000 else f"t={t / 1440:8.1f} d"
        print(f"  {label}  " + "  ".join(parts))
    if args.depth > 0:
        for line in body.tree(depth=args.depth):
            print("  " + line)
    if body.spatial:
        print(f"  space {o.width}x{o.height}; blocked divisions (no free site): {body.blocked_divisions}")
        for row in body.type_map():
            print("  " + row)
        bands = body.bands_along_x()
        if bands:
            print("  bands: " + ", ".join(f"{t or 'empty'} {a}-{b}" for t, a, b in bands))
    if body.unknown:
        print("  UNKNOWN stops: " + ", ".join(f"{k} ×{v}" for k, v in body.unknown.items()))
    for chk in body.check_asserts():
        mark = "ok  " if chk["ok"] else "FAIL"
        print(f"  assert {mark} {chk['assert']}  (got {chk.get('value', chk.get('error'))})")
    print(body.uncertainty().format())
    out: dict = {"summary": s, "asserts": body.check_asserts(), "uncertainty": body.uncertainty().to_dict()}
    if args.compare:
        path = REFERENCES.get(o.reference)
        if not path or not Path(path).exists():
            print(f"  no reference {o.reference!r} available (genomeos data distil --only celegans_lineage)")
        else:
            diff = compare(body, ReferenceLineage.load(path), until=args.until)
            print(diff.format())
            out["diff"] = diff.to_dict()
    if args.json:
        Path(args.json).write_text(json.dumps(out, indent=1, default=str))
        print(f"  wrote {args.json}")
    return 0


def cmd_data(args: argparse.Namespace) -> int:
    from genomeos import storage
    from genomeos.results import list_results

    if args.data_cmd == "fetch":
        from genomeos.genome.fetch import fetch_chromosome, fetch_individual, individual_vcf_path

        if args.individual:
            missing = [c for c in args.chrom if not individual_vcf_path(c).exists()]
            if missing:
                print("streaming the GIAB HG002 benchmark once (156 MB) and keeping every chromosome's rows…")
                counts = fetch_individual(progress=lambda m: print(f"  {m}", flush=True))
                print(f"  kept {sum(counts.values()):,} PASS variants over {len(counts)} chromosomes")
            for c in args.chrom:
                print(f"  {c}: {individual_vcf_path(c)}")

        for chrom in args.chrom:
            print(f"{chrom}: sequence, gene models, regulatory elements, repeats…", flush=True)
            r = fetch_chromosome(
                chrom,
                progress=lambda m: print(f"  {m}", flush=True),
                elements=not args.no_elements,
                repeats=not args.no_repeats,
                analyse=args.analyse,
            )
            print(
                f"  ready: {r['sequence']}, {r['gencode']}"
                + (f", {r['ccres']} ENCODE elements" if "ccres" in r else "")
                + (f", {r['repeats']:,} repeat copies" if "repeats" in r else "")
            )
        print("every view now works for these chromosomes: Blocks, Flow, regulation, domains, unknown")
        return 0
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


def cmd_calibrate(args: argparse.Namespace) -> int:
    from genomeos.twin.calibrate import calibrate_attrition

    cal = calibrate_attrition(args.cell_type, args.age, args.measured, args.source)
    print(f"{args.cell_type}: measured {args.measured:.0f} bp at {args.age:.0f} y")
    print(f"  fitted divisions_per_year = {cal.fitted_divisions_per_year:.3f}")
    print(f"  net attrition {cal.net_bp_per_year:.1f} bp/yr")
    print(f"  saved as calibration_{args.cell_type}_attrition; twin runs of this cell type now use it")
    ev = cal.parameter.evidence
    print(f"  [{ev.kind.value}] {ev.source}  conf={cal.parameter.confidence}")
    return 0


def cmd_methylation(args: argparse.Namespace) -> int:
    from genomeos.twin.clocks import Clock
    from genomeos.twin.methylation import load_probes, summarise

    probes = load_probes()
    summary = summarise(
        args.sources,
        probes,
        min_coverage=args.min_coverage,
        progress=lambda n: print(f"  {n:,} rows scanned", flush=True),
    )
    print(f"{summary.rows_scanned:,} bedMethyl rows scanned")
    print(
        f"  {summary.probes_covered}/{summary.probes_total} clock CpGs covered (>= {args.min_coverage} reads)"
    )
    for clock in (Clock.horvath(), Clock.hannum()):
        r = clock.predict(summary.betas)
        conf = r.parameter.confidence
        print(f"  {clock.name:<12} {r.age:5.1f} y  (coverage {r.coverage:.0%}, conf {conf:.2f})")
    if args.save:
        from genomeos.results import save_result

        h = Clock.horvath().predict(summary.betas)
        p = save_result(
            args.save,
            {
                "sources": args.sources,
                "rows_scanned": summary.rows_scanned,
                "probes_covered": summary.probes_covered,
                "horvath_age": round(h.age, 2),
                "hannum_age": round(Clock.hannum().predict(summary.betas).age, 2),
                "betas": {k: round(v, 4) for k, v in summary.betas.items()},
                "disk_used_bytes": 0,
            },
        )
        print(f"  saved {p}")
    return 0


def cmd_signals(args: argparse.Namespace) -> int:
    from genomeos.genome import Annotation, IndexedGenome, default_gencode, learn_signals, scan
    from genomeos.genome.signals import SignalSet

    if args.signals_cmd == "learn":
        gff = args.gff3 or default_gencode({args.chrom})
        ann = Annotation.from_gff3(gff, {args.chrom})
        genome = IndexedGenome(args.genome)
        sig = learn_signals(ann, genome, args.chrom)
        genome.close()
        out = Path(args.output or f"data/results/signals_{args.chrom}.json")
        sig.save(out)
        st = sig.stats
        print(f"learned from {st['transcripts']} transcripts, {st['introns']} introns on {args.chrom}")
        for k, v in sig.stats.items():
            if k not in ("chromosome", "transcripts", "introns"):
                print(f"  {k:<28} {v}")
        for name, p in sig.pwms.items():
            print(f"  {name:<18} consensus {p.consensus}  ({p.examples} examples)")
        print(f"  saved {out}")
        return 0
    sig = SignalSet.load(args.signals)
    genome = IndexedGenome(args.genome)
    from genomeos.coords import Locus

    loc = Locus.parse(args.locus)
    seq = genome.fetch(Locus(loc.chrom, loc.start, loc.end))
    genome.close()
    hits = scan(seq, sig, min_relative=args.min_relative)
    print(f"{loc}: {len(hits)} signal hits at relative score >= {args.min_relative}")
    for h in hits[: args.limit]:
        where = f"{loc.chrom}:{loc.start + h.pos}"
        print(f"  {where:<20} {h.strand.value}  {h.signal:<16} score {h.score:6.2f}  rel {h.relative:.2f}")
    return 0


def _anatomy_table(anatomies: list) -> str:
    rows = []
    keys = [
        ("length (Mb)", lambda a: f"{a.length / 1e6:.3f}"),
        ("coding genes", lambda a: a.genes.get("protein_coding", "-")),
        ("genes per Mb", lambda a: a.layout.get("coding_genes_per_mb", "-")),
        ("CDS fraction", lambda a: f"{a.composition_fraction().get('cds', 0):.1%}"),
        ("intron fraction", lambda a: f"{a.composition_fraction().get('intron', 0):.1%}"),
        ("intergenic fraction", lambda a: f"{a.composition_fraction().get('intergenic', 0):.1%}"),
        ("exons / transcript", lambda a: a.structure.get("exons_per_transcript_median", "-")),
        ("intron median (bp)", lambda a: a.structure.get("intron_length_median", "-")),
        ("gene median (bp)", lambda a: a.structure.get("gene_length_median", "-")),
        ("protein median (aa)", lambda a: a.structure.get("protein_length_aa_median", "-")),
        ("single-exon genes", lambda a: a.structure.get("single_exon_fraction", "-")),
        ("transcripts / gene", lambda a: a.genes.get("transcripts_per_coding_gene_mean", "-")),
        ("spacing median (bp)", lambda a: a.layout.get("intergenic_spacing_median", "-")),
        ("overlapping pairs", lambda a: a.layout.get("overlapping_coding_pairs", "-")),
        ("GC", lambda a: a.sequence.get("gc", "-")),
        ("CpG islands", lambda a: a.elements.get("cpg_islands", "-")),
        ("homopolymers >=12", lambda a: a.elements.get("homopolymer_runs_ge12", "-")),
    ]
    for label, fn in keys:
        rows.append({"measure": label, **{a.name: fn(a) for a in anatomies}})
    return _table(rows, ["measure", *[a.name for a in anatomies]])


def cmd_anatomy(args: argparse.Namespace) -> int:
    from genomeos.genome import Annotation, Genome, anatomy_of, default_gencode, design_lessons
    from genomeos.results import load_result, save_result

    if args.compare:
        from genomeos.genome.anatomy import Anatomy

        r = load_result("anatomy_comparison")
        if not r:
            print("no saved comparison; run `genomeos anatomy` on genomes with --save first")
            return 1
        ans = [Anatomy(**{k: v for k, v in d.items() if k != "composition_fraction"}) for d in r["anatomies"]]
        print(_anatomy_table(ans))
        print()
        for line in r["lessons"]:
            print("  " + line)
        return 0
    genome = Genome.from_fasta(args.fasta)
    chrom = args.chrom or next(iter(genome.chromosomes))
    seq = genome.chromosomes[chrom].sequence
    ann = None
    gff = args.gff3 or (default_gencode({chrom}) if chrom.startswith("chr") else None)
    if gff:
        ann = Annotation.from_gff3(gff, {chrom})
    a = anatomy_of(args.name or f"{Path(args.fasta).stem} {chrom}", seq, ann, chrom if ann else None)
    d = a.to_dict()
    print(f"{a.name}: {a.length:,} bp")
    print("composition:")
    for k, v in a.composition_bp.items():
        bar = "█" * int(v / a.length * 40)
        print(f"  {k:<24} {v:>12,} bp  {v / a.length:6.1%}  {bar}")
    for section in ("sequence", "elements", "genes", "structure", "layout"):
        if d[section]:
            print(f"{section}:")
            for k, v in d[section].items():
                print(f"  {k:<36} {v}")
    if ann:
        for line in design_lessons([a]):
            print("lesson: " + line)
    if args.save:
        save_result(args.save, d)
        print(f"saved data/results/{args.save}.json")
    return 0


def cmd_design(args: argparse.Namespace) -> int:
    from genomeos.design import IDENTITY, design

    if args.cell_type not in IDENTITY:
        print(f"known cell types: {', '.join(sorted(IDENTITY))}")
        return 1
    d = design(args.cell_type)
    print(f"design budget for a {args.cell_type}")
    print(f"  essential genes (every cell, CEGv2)       {d.essential_genes:>7,}")
    print(f"  identity genes (chosen libraries − core) {d.identity_genes:>7,}")
    masters = ", ".join(d.master_regulators)
    print(f"  master regulators                         {len(d.master_regulators):>7}   {masters}")
    print(f"  minimal program (essential+identity+masters) {d.minimal_genes:>4,} genes")
    print(f"  broad program (all core libraries+identity) {d.total_genes:>5,} genes")
    print(f"    ({d.span_bp / 1e6:.0f} Mb of human gene spans)")
    print("  base budget for the minimal program:")
    for k, v in d.budgets_bp.items():
        print(f"    {k:<40} {v / 1e6:8.1f} Mb")
    print("  libraries used:")
    for k, v in d.libraries.items():
        print(f"    {k:<40} {v:>6,} genes")
    for n in d.notes:
        print(f"  note: {n}")
    if args.save:
        from genomeos.results import save_result

        print(f"  saved {save_result(f'design_{args.cell_type}', d.to_dict())}")
    return 0


def cmd_cancer(args: argparse.Namespace) -> int:
    from genomeos.results import load_result

    k = load_result("cancer_msk_impact_2017")
    if args.cancer_cmd == "distil":
        from genomeos.cancer import CBioPortal
        from genomeos.results import save_result

        d = CBioPortal().distil_study(args.study, progress=lambda i, n: print(f"  {i}/{n} genes", flush=True))
        print(f"  saved {save_result(f'cancer_{args.study}', d)}")
        return 0
    if args.cancer_cmd == "expression":
        from genomeos.cancer import CBioPortal
        from genomeos.cancer.cbioportal import DRIVER_PANEL
        from genomeos.results import save_result

        genes = args.genes or list(DRIVER_PANEL)
        d = CBioPortal(timeout=180).expression_distribution(
            args.study, genes, progress=lambda i, n: print(f"  {i}/{n} genes", flush=True)
        )
        print(f"{d['study']}: {d['samples']} tumours, {len(d['genes'])} genes, profile {d['kind']}")
        print(f"  {d['meaning']}")
        rows = sorted(d["genes"].items(), key=lambda kv: -(kv[1].get("fraction_raised") or kv[1]["median"]))[
            : args.top
        ]
        print(
            _table(
                [
                    {
                        "gene": g,
                        "median": f"{v['median']:+.2f}",
                        "q1..q3": f"{v['q1']:+.2f}..{v['q3']:+.2f}",
                        "raised above normal": (
                            f"{v['fraction_raised']:.1%}" if "fraction_raised" in v else "-"
                        ),
                    }
                    for g, v in rows
                ],
                ["gene", "median", "q1..q3", "raised above normal"],
            )
        )
        if d["missing"]:
            print(f"  not measured: {', '.join(d['missing'][:12])}")
        print(f"  saved {save_result(f'expression_{args.study}', d)}")
        return 0
    if not k:
        print("no distilled cancer knowledge yet: run `genomeos cancer distil` (cBioPortal, ~30 s)")
        return 1
    if args.cancer_cmd == "genes":
        rows = sorted(k["genes"].items(), key=lambda kv: -kv[1]["frequency"])[: args.top]
        print(f"{k['study']}: {k['samples']:,} tumours, {len(k['genes'])} driver genes")
        print(
            _table(
                [
                    {
                        "gene": g,
                        "mutated": f"{v['frequency']:.1%}",
                        "top change": v["hotspots"][0][0] if v["hotspots"] else "",
                        "types (>10%)": ", ".join(
                            f"{t} {f:.0%}"
                            for t, f in sorted(v["by_cancer_type"].items(), key=lambda x: -x[1])
                            if f > 0.1
                        )[:70],
                    }
                    for g, v in rows
                ],
                ["gene", "mutated", "top change", "types (>10%)"],
            )
        )
        return 0
    if args.cancer_cmd == "gene":
        v = k["genes"].get(args.symbol)
        if not v:
            print(f"{args.symbol} not in the distilled panel")
            return 1
        print(f"{args.symbol}: mutated in {v['frequency']:.1%} of {k['samples']:,} tumours ({k['study']})")
        print("  hotspots:", ", ".join(f"{h[0]} ×{h[1]}" for h in v["hotspots"]))
        for ct, f in sorted(v["by_cancer_type"].items(), key=lambda x: -x[1])[:12]:
            print(f"  {ct:<36} {f:6.1%}")
        return 0
    if args.cancer_cmd == "types":
        print(
            _table(
                [{"cancer type": t, "samples": n} for t, n in list(k["cancer_types"].items())[: args.top]],
                ["cancer type", "samples"],
            )
        )
        return 0
    if args.cancer_cmd == "tumour":
        from genomeos.cancer import analyse, tumour_packet

        a = analyse(
            args.vcf,
            set(args.chrom) if args.chrom else None,
            k,
            deep=args.deep,
            pathways=not args.no_pathways,
            therapeutics=args.therapeutic,
            log=sys.stdout,
        )
        mb = a["mutation_burden"]
        print(
            f"{args.vcf}: {a['variants_total']} variants, {a['annotated']} annotated by VEP; "
            f"{a['likely_germline']} set aside as likely germline; {a['somatic_candidates']} somatic "
            f"candidates, {a['coding_somatic']} coding"
        )
        print(f"  mutation burden: {mb['per_megabase']} coding/Mb [{mb['evidence']}; {mb['assumption']}]")
        print(
            _table(
                [
                    {
                        "gene": r["gene"] or "-",
                        "consequence": r["consequence"],
                        "change": r["protein_change"] or "-",
                        "driver": f"{r['driver_frequency']:.1%}" if r["driver_frequency"] else "-",
                        "hotspot": "yes" if r["hotspot"] else "",
                        "COSMIC": "yes" if r["cosmic"] else "",
                        "damage": r["sift"] or r["polyphen"] or "-",
                        "score": f"{r['score']:.2f}",
                    }
                    for r in a["ranked"][:15]
                ],
                ["gene", "consequence", "change", "driver", "hotspot", "COSMIC", "damage", "score"],
            )
        )
        if a["germline_set_aside"]:
            print(
                "  likely germline (population frequency): "
                + ", ".join(
                    f"{g['gene'] or g['chrom'] + ':' + str(g['pos'])} {g['gnomad_af']:.2g}"
                    for g in a["germline_set_aside"][:8]
                )
            )
        for p in a["mutant_peptides"]:
            if "error" in p:
                print(f"  peptide {p['error']}")
            else:
                print(
                    f"  peptide {p['gene']} {p['change']}: {p['wild_type']} → {p['mutant']} "
                    f"(residues {p['window'][0]}-{p['window'][1]})"
                )
        for e in a["pathway_effects"]:
            if "error" in e:
                print(f"  pathway {e['gene']}: {e['error']}")
            else:
                w = e["most_affected"]
                print(
                    f"  pathway {e['gene']} {e['change']} treated as absent: {e['reactions_lost']} reactions "
                    f"lost over {e['pathways_checked']} pathways"
                    + (f"; most affected {w['name']} ({w['fraction_lost']:.0%})" if w else "")
                )
        print("suggested cancer types (inferred):")
        for t_ in a["suggested_cancer_types"]:
            print(f"  {t_['cancer_type']:<36} enrichment {t_['log_enrichment']:+.2f}")
        print(
            "cell-surface products among altered genes:",
            ", ".join(x["gene"] for x in a["surface_targets"]) or "none",
        )
        t = a.get("therapeutic_targets")
        if t:
            print(f"therapeutic target candidates (data level {t['data_level']}/9):")
            for c in t["therapeutic_candidates"]:
                best = (c["recommended_mechanisms"] or [{}])[0]
                print(
                    f"  {c['gene']:<10} {c['target_class']:<24} score "
                    f"{c['scores']['overall'] if c['scores']['overall'] is not None else '-'}"
                    + (
                        f"  best mechanism {best['mechanism']} "
                        f"{best['compatibility']:.0%} (cargo {best['cargo']})"
                        if best
                        else "  no mechanism passes its requirements"
                    )
                )
            print(f"  {t['disclaimer']}")
        if args.packet:
            Path(args.packet).write_text(json.dumps(tumour_packet(a), indent=1))
            print(f"agent packet written to {args.packet}")
        return 0
    # compare
    from genomeos.cancer import agent_packet, annotate, somatic, suggest_cancer_type, surface_targets
    from genomeos.genome import Annotation, IndexedGenome, default_gencode

    chroms = set(args.chrom) if args.chrom else None
    som = somatic(args.normal, args.tumour, chroms)
    print(f"somatic variants (tumour − normal): {len(som)}")
    ann = Annotation.from_gff3(args.gff3 or default_gencode(chroms), chroms)
    genome = IndexedGenome(args.genome)
    ranked = annotate(som, ann, genome, k)
    genome.close()
    print(
        _table(
            [
                {
                    "gene": s.gene or "-",
                    "consequence": s.consequence or "-",
                    "change": s.protein_change,
                    "driver freq": f"{s.driver_frequency:.1%}" if s.driver_frequency else "-",
                    "hotspot": "yes" if s.hotspot else "",
                    "score": f"{s.score:.2f}",
                }
                for s in ranked[:15]
            ],
            ["gene", "consequence", "change", "driver freq", "hotspot", "score"],
        )
    )
    genes = {
        s.gene
        for s in ranked
        if s.gene and s.consequence not in ("synonymous_variant", "intron_variant", "intergenic")
    }
    types = suggest_cancer_type(genes, k)
    print("suggested cancer types (inferred):")
    for t_ in types:
        print(f"  {t_['cancer_type']:<36} enrichment {t_['log_enrichment']:+.2f}")
    targets = surface_targets(genes)
    print("cell-surface products among altered genes:", ", ".join(x["gene"] for x in targets) or "none")
    if args.packet:
        Path(args.packet).write_text(json.dumps(agent_packet(args.tumour, ranked, types, targets), indent=1))
        print(f"agent packet written to {args.packet}")
    return 0


def _pct(x: float | None, digits: int = 1) -> str:
    return f"{x:.{digits}%}" if x is not None else ""


def cmd_budget(args: argparse.Namespace) -> int:
    """The 98%: every UNKNOWN block with constraint attached and a best guess (docs/ATTRIBUTION.md)."""
    import time

    from genomeos.attribution.budget import TIERS, distil, run_and_save
    from genomeos.results import load_result, save_result

    if args.chrom and args.bio:
        from genomeos.attribution.compile import compile_chromosome, write_program

        if args.out:
            p = write_program(args.chrom, Path(args.out))
            print(f"wrote {p} ({p.stat().st_size / 1e3:.0f} kB); check it with: bio test {p}")
        else:
            print(compile_chromosome(args.chrom), end="")
        return 0
    if args.chrom:
        out = load_result(f"budget_{args.chrom}") if args.cached else None
        if out is None:
            last = [0.0]

            def progress(done: int, total: int, mb: int) -> None:
                if time.time() - last[0] > 15:
                    last[0] = time.time()
                    print(f"  phyloP: {done}/{total} data blocks, {mb / 1e6:.0f} MB", flush=True)

            out = run_and_save(
                args.chrom,
                threshold=args.threshold,
                phylop=not args.no_phylop,
                elements=not args.no_elements,
                progress=progress,
            )
        cf = out.get("constrained_fraction")
        print(
            f"{args.chrom}: {len(out['blocks'])} UNKNOWN blocks, {out['unknown_bp'] / 1e6:.1f} Mb of "
            f"{out['chromosome_length'] / 1e6:.1f} Mb; constrained bases (phyloP >= {out['threshold']}) "
            + (f"{cf:.2%} of those measured" if cf is not None else "not measured")
        )
        rows = [
            {
                "tier": t,
                "blocks": v["blocks"],
                "Mb": f"{v['bp'] / 1e6:.2f}",
                "of UNKNOWN": f"{v['fraction_of_unknown']:.1%}",
                "of chromosome": _pct(v.get("fraction_of_chromosome")),
                "constrained kb": f"{v['constrained_bp'] / 1e3:.1f}",
            }
            for t, v in out["by_tier"].items()
        ]
        print(_table(rows, ["tier", "blocks", "Mb", "of UNKNOWN", "of chromosome", "constrained kb"]))
        top = sorted(
            (r for r in out["blocks"] if r["guess"]["tier"] == "constrained_unknown"),
            key=lambda r: -((r["phylop"] or {}).get("above") or 0),
        )[: args.top]
        if top:
            print(f"\nconstrained_unknown, most constrained bases first (top {len(top)}):")
            print(
                _table(
                    [
                        {
                            "block": f"{r['start']:,}-{r['end']:,}",
                            "kb": f"{r['length'] / 1e3:.0f}",
                            "class": r["class"],
                            "constrained": f"{(r['phylop'] or {}).get('fraction_above', 0) or 0:.1%}",
                            "elements": (r["elements"] or {}).get("n", ""),
                            "max lod": (r["elements"] or {}).get("max_lod", ""),
                        }
                        for r in top
                    ],
                    ["block", "kb", "class", "constrained", "elements", "max lod"],
                )
            )
        print(
            "  [phyloP: Zoonomia 241 mammals, read per base from UCSC, never stored; elements: 100"
            " vertebrates; a tier is a best guess with its confidence, not a verdict]"
        )
        return 0
    s = distil()
    if not s["chromosomes"]:
        print("no budget_chr*.json yet: genomeos budget --chrom chr21, or the budget_genome_wide job")
        return 1
    save_result("budget_genome_wide", s)
    print(
        f"{s['chromosomes']} chromosomes: {s['unknown_bp'] / 1e6:.0f} Mb of UNKNOWN blocks in "
        f"{s['genome_bp'] / 1e6:.0f} Mb; constrained {s['constrained_fraction']:.2%} of measured bases; "
        f"a guess at confidence >= 0.5 on {s['guessed_fraction']:.1%}"
    )
    rows = [
        {
            "tier": t,
            "blocks": v["blocks"],
            "Mb": f"{v['bp'] / 1e6:.1f}",
            "of UNKNOWN": _pct(v.get("fraction_of_unknown")),
            "of genome": _pct(v.get("fraction_of_genome")),
            "constrained Mb": f"{v['constrained_bp'] / 1e6:.2f}",
        }
        for t, v in s["by_tier"].items()
        if t in TIERS
    ]
    print(_table(rows, ["tier", "blocks", "Mb", "of UNKNOWN", "of genome", "constrained Mb"]))
    rows = [
        {
            "class": k,
            "blocks": v["blocks"],
            "Mb": f"{v['bp'] / 1e6:.1f}",
            "constrained": _pct(v.get("constrained_fraction"), 2),
        }
        for k, v in s["by_class"].items()
    ][: args.top]
    print(_table(rows, ["class", "blocks", "Mb", "constrained"]))
    return 0


def cmd_closure(args: argparse.Namespace) -> int:
    """The gene-level closure test: do the attributions reproduce what each cell makes? (area I)"""
    import time

    from genomeos.attribution.closure import CELLS, run_and_save
    from genomeos.results import load_result

    out = load_result(f"closure_{args.chrom}") if args.cached else None
    if out is None:
        last = [0.0]

        def progress(done: int, total: int, mb: int) -> None:
            if time.time() - last[0] > 15:
                last[0] = time.time()
                print(f"  RNA: {done}/{total} data blocks, {mb / 1e6:.0f} MB", flush=True)

        cells = tuple(c.strip() for c in args.cells.split(",")) if args.cells else CELLS
        out = run_and_save(args.chrom, cells=cells, signal=args.signal, progress=progress)
    swapped = [k for k, v in out.get("strand_orientation", {}).items() if v == "swapped"]
    print(
        f"{args.chrom}: {out['coding_genes']} coding genes, {out['genes_with_elements']} with attributed "
        f"elements, cells {', '.join(out['cells'])}; expressed = canonical exons "
        f"≥ {out['expressed_threshold']:.0%} covered at signal {out['signal_threshold']}"
        + (f"; RNA strands read swapped for {', '.join(swapped)}" if swapped else "")
    )

    def f(w: dict, k: str) -> str:
        return f"{w[k][0]:.0%} ({w[k][1]})" if w[k][0] is not None else f"- ({w[k][1]})"

    rows = []
    for cell, w in out["within_cell"].items():
        rows.append(
            {
                "cell": cell,
                "promoter open": f(w, "expressed_promoter_open"),
                "promoter closed": f(w, "expressed_promoter_closed"),
                "open + activating input": f(w, "expressed_open_with_activating_input"),
                "open, no element": f(w, "expressed_open_without_input"),
                "open + repressing input": f(w, "expressed_open_with_repressing_input"),
            }
        )
    print("\nexpressed fraction (genes) by what the attribution says:")
    print(
        _table(
            rows,
            [
                "cell",
                "promoter open",
                "promoter closed",
                "open + activating input",
                "open, no element",
                "open + repressing input",
            ],
        )
    )
    a = out["across_cells"]
    print(
        f"\nacross cells, {a['genes_tested']} genes with elements and variation: most active cell is "
        f"the most expressed for {_pct(a['most_active_cell_is_most_expressed'], 0)} (shuffled cells "
        f"{_pct(a['same_under_shuffled_cells'], 0)}, chance {_pct(a['chance'], 0)}); mean rank correlation "
        f"input vs expression {a['mean_rho_input_vs_expression']}, promoter vs expression "
        f"{a['mean_rho_promoter_vs_expression']}"
    )
    for mode, m in (out.get("modes") or {}).items():
        if mode == "dnase":
            continue
        a = m["across_cells"]
        label = {
            "cell": "the deletion scored on the cell's own track, no DNase gating",
            "both": "the cell's own score, only where a DNase peak sits on the element",
        }.get(mode, mode)
        print(
            f"  with {label}: {a['genes_tested']} genes, most active cell is the most expressed for "
            f"{_pct(a['most_active_cell_is_most_expressed'], 0)} (shuffled "
            f"{_pct(a['same_under_shuffled_cells'], 0)}), mean rank correlation "
            f"{a['mean_rho_input_vs_expression']}; rejected {m['rejected_count']}"
        )
    if out["rejected_attributions"]:
        print(
            f"\nrejected by the cell ({out['rejected_count']}): promoter open, activating input, gene silent"
        )
        print(_table(out["rejected_attributions"][: args.top], ["gene", "cell", "input", "expression"]))
    if out["unexplained_expression"]:
        print(f"\nunexplained ({out['unexplained_count']}): expressed, promoter closed, no active element")
        print(_table(out["unexplained_expression"][: args.top], ["gene", "cell", "expression"]))
    print(
        "  [expression: ENCODE total RNA-seq over canonical exons; promoter: DNase peak at the TSS; "
        "elements: AlphaGenome deletion target, active where a DNase peak overlaps; docs/ATTRIBUTION.md]"
    )
    return 0


def cmd_decompile(args: argparse.Namespace) -> int:
    """One gene, every layer GenomeOS has read for it, as a BioLang-flavoured view (area J step 6)."""
    from genomeos.decompile import decompile, render

    d = decompile(args.gene, args.chrom)
    if args.json:
        print(json.dumps(d, indent=1))
        return 0
    print(render(d), end="")
    return 0


def cmd_motifs(args: argparse.Namespace) -> int:
    """JASPAR motifs over promoters and elements: `requires:` lists and operators (area J step 4)."""
    from genomeos.results import load_result

    if args.library:
        r = load_result("motifs_genome_wide")
        if r is None:
            print("no motifs_genome_wide result yet; run: uv run python scripts/motifs_genome_wide.py")
            return 1
        lib = (r.get("libraries") or {}).get(args.library)
        if not lib:
            print(f"{args.library}: fewer members with a scanned promoter than the summary needs")
            return 1
        print(f"{args.library}: {lib['members_placed']} of {lib['members']} members with a scanned promoter")
        if lib["factors"]:
            print(
                _table(
                    lib["factors"][: args.top], ["factor", "members_with_hit", "share", "all_share", "ratio"]
                )
            )
        if lib["operators"]:
            print("\noperators (factor pairs recurring across the members beyond the two shares):")
            print(
                _table(
                    [{**o, "factors": " + ".join(o["factors"])} for o in lib["operators"][: args.top]],
                    ["factors", "members_with_both", "share", "expected", "ratio"],
                )
            )
        print(f"  [{r['evidence']['operator']}]")
        return 0
    if not args.chrom:
        print("give --chrom C (a chromosome's promoters and elements) or --library L (the genome summary)")
        return 1
    r = load_result(f"motifs_{args.chrom}") if args.cached else None
    if r is None:
        from genomeos.genome.motifs import run_and_save

        r = run_and_save(args.chrom, progress=lambda m: print(f"  {m}", flush=True))
    if args.gene:
        g = r["genes"].get(args.gene)
        if not g:
            print(f"{args.gene}: no scanned promoter on {args.chrom}")
            return 1
        print(
            f"{args.gene}: {g['factors_hit']} factors hit the promoter; "
            "requires (enriched, best score first):"
        )
        print(_table(g["requires"], ["factor", "score", "position", "strand", "enrichment"]))
        print(f"  [{r['evidence']['requires']}]")
        return 0
    print(
        f"{args.chrom}: {r['promoters']} promoters (TSS ± {r['flank']} bp) and {len(r['elements'])} elements "
        f"scanned against {r['profiles']} JASPAR profiles ({r['factors']} factors); "
        f"{r['requires_per_promoter']} enriched factors per promoter on average"
    )
    print(_table(r["most_enriched"][: args.top], ["factor", "sequences", "shuffled", "share", "enrichment"]))
    libs = [(k, v) for k, v in r["libraries"].items() if v["factors"]]
    if libs:
        print("\nlibraries whose members' promoters share a factor beyond the chromosome's rate:")
        for k, v in libs[: args.top]:
            print(
                f"  {k} ({v['members_on_chromosome']} members): "
                + ", ".join(f"{f['factor']} {f['ratio']}x" for f in v["factors"][:5])
            )
    print(f"  [{r['evidence']['hit']}; {r['cost']['seconds']} s]")
    return 0


def cmd_duplications(args: argparse.Namespace) -> int:
    """Segmental duplications over the UNKNOWN blocks: curated copy-and-paste (area J step 3)."""
    from genomeos.genome.duplications import run_and_save
    from genomeos.results import load_result

    r = load_result(f"duplication_{args.chrom}") if args.cached else None
    if r is None:
        r = run_and_save(args.chrom)
    print(
        f"{args.chrom}: {r['pairs']:,} curated duplication pairs ({r['intra_chromosomal']:,} within the "
        f"chromosome), {r['duplicated_bp'] / 1e6:.2f} Mb duplicated "
        + (
            f"({_pct(r['duplicated_fraction'])} of the chromosome)"
            if r["duplicated_fraction"] is not None
            else ""
        )
        + f"; identity bands young {r['bands']['young']}, middle {r['bands']['middle']}, "
        f"old {r['bands']['old']}"
    )
    if r["by_tier"]:
        print(
            _table(
                [
                    {
                        "tier": t,
                        "blocks": v["blocks"],
                        "Mb": f"{v['bp'] / 1e6:.2f}",
                        "duplicated": _pct(v["duplicated_fraction"]),
                        "mostly duplicated": v["blocks_mostly_duplicated"],
                    }
                    for t, v in r["by_tier"].items()
                ],
                ["tier", "blocks", "Mb", "duplicated", "mostly duplicated"],
            )
        )
    sp = r.get("similar_to_pairs")
    if sp:
        print(
            f"  the classifier's shared-k-mer pairs: {sp['pairs']} pairs, "
            f"{sp['supported_by_curated_duplication']} supported by a curated duplication "
            f"({_pct(sp['share'])})"
        )
    el = r.get("elements")
    if el:
        print(
            f"  scored elements inside a duplication: {el['inside_a_duplication']} of {el['scored']} "
            f"({_pct(el['share'])})"
        )
    top = sorted((b for b in r["blocks"] if b["pairs"]), key=lambda b: -b["duplicated_fraction"])[: args.top]
    if top:
        print(f"\nmost duplicated UNKNOWN blocks (top {len(top)}):")
        print(
            _table(
                [
                    {
                        "locus": f"{args.chrom}:{b['start']:,}-{b['end']:,}",
                        "kb": f"{(b['end'] - b['start']) / 1e3:.0f}",
                        "tier": b["tier"],
                        "duplicated": _pct(b["duplicated_fraction"]),
                        "pairs": b["pairs"],
                        "partner": b["partners"][0] if b["partners"] else "",
                    }
                    for b in top
                ],
                ["locus", "kb", "tier", "duplicated", "pairs", "partner"],
            )
        )
    print(f"  [{r['evidence']['pairs']}; {r['cost']['seconds']} s]")
    return 0


def cmd_across(args: argparse.Namespace) -> int:
    """One locus across species: alignment, conserved grammar, both constraint axes (area J step 7)."""
    from genomeos.knowledge.across import LOCI
    from genomeos.results import load_result

    if args.locus not in LOCI:
        print(f"unknown locus {args.locus}; known: {', '.join(LOCI)}")
        return 1
    r = load_result(f"across_{args.locus}") if args.cached else None
    if r is None:
        from genomeos.knowledge.across import run_and_save

        r = run_and_save(args.locus, progress=lambda m: print(f"  {m}", flush=True))
    print(
        f"{args.locus}: {r['chrom']}:{r['start']:,}-{r['end']:,} ({r['length']} bp), "
        f"target {r['gene']}; {r['note']}"
    )
    print(
        _table(
            [
                {
                    "species": sp,
                    "aligned": _pct(v.get("coverage")),
                    "identity": _pct(v.get("identity")),
                    "locus": v.get("locus") or v.get("error") or "no alignment",
                }
                for sp, v in r["species"].items()
            ],
            ["species", "aligned", "identity", "locus"],
        )
    )
    c = r["constraint"]
    ph, gn = c["phylop"], c["gnocchi"]
    print(
        f"  mammals: phyloP mean {ph['mean']}, {_pct(ph['fraction_above'])} of bases constrained; "
        + (
            f"people: {_pct(gn['fraction_above'])} of kilobases constrained (max Z {gn['maximum']})"
            if gn
            else "people: kilobase not scored by Gnocchi"
        )
    )
    g = r.get("grammar") or {}
    if g:
        print(
            f"  {g['factors_in_human']} factors hit the human element; holding the same aligned site in "
            f"every species: {', '.join(g['everywhere']) or 'none'}"
        )
        rows = [
            {
                "factor": x["factor"],
                "human": x["human_score"],
                "same site in": ", ".join(x["same_site_in"]) or "human only",
            }
            for x in g["rows"][: args.top]
        ]
        print(_table(rows, ["factor", "human", "same site in"]))
    print(f"  [{r['evidence']['alignment']}; {r['evidence']['grammar']}; {r['cost']['seconds']} s]")
    return 0


def cmd_profiling(args: argparse.Namespace) -> int:
    """Phylogenetic profiling between libraries: gained and lost together across species (area J step 5)."""
    from genomeos.knowledge.homology import STRATA
    from genomeos.results import load_result

    r = load_result("profiling_genome_wide") if args.cached else None
    if r is None:
        from genomeos.knowledge.profiling import run_and_save

        r = run_and_save()
    if args.library:
        lib = r["libraries"].get(args.library)
        if not lib:
            print(f"{args.library}: fewer than {r['min_members']} members placed")
            return 1
        print(
            f"{args.library}: {lib['members_placed']} members; profile against the genome's r = "
            f"{lib['against_genome']}; gain curve (share of members with an origin at or before the stratum):"
        )
        gc = lib["gain_curve"]
        print(
            _table(
                [{"stratum": s_, "share": _pct(gc[s_])} for s_ in STRATA if s_ in gc], ["stratum", "share"]
            )
        )
        partners = [p_ for p_ in r["pairs"] if args.library in p_["libraries"]]
        rows = [
            {"library": [x for x in p_["libraries"] if x != args.library][0], "r": p_["residual_correlation"]}
            for p_ in partners[: args.top]
        ]
        if rows:
            print("\nlibraries with the most similar history (residual profile correlation):")
            print(_table(rows, ["library", "r"]))
        print(f"  [{r['evidence']['correlation']}]")
        return 0
    print(
        f"{len(r['libraries'])} libraries profiled over {r['species']} species "
        f"({', '.join(f'{k} {v}' for k, v in r['species_by_stratum'].items())})"
    )
    print(
        _table(
            [
                {"libraries": " ~ ".join(p_["libraries"]), "r": p_["residual_correlation"]}
                for p_ in r["top_pairs"][: args.top]
            ],
            ["libraries", "r"],
        )
    )
    print("\nleast alike:")
    print(
        _table(
            [
                {"libraries": " ~ ".join(p_["libraries"]), "r": p_["residual_correlation"]}
                for p_ in r["bottom_pairs"][-min(args.top, 5) :]
            ],
            ["libraries", "r"],
        )
    )
    print(f"  [{r['evidence']['profile']}; {r['evidence']['correlation']}]")
    return 0


def cmd_origin(args: argparse.Namespace) -> int:
    """Origin per gene and age per library from Ensembl Compara, streamed once (area J step 2)."""
    from genomeos.knowledge.homology import LADDER, STRATA
    from genomeos.results import load_result

    r = load_result("origin_genome_wide")
    if r is None:
        print("no origin_genome_wide result yet; run: uv run python scripts/origin_genome_wide.py")
        return 1
    if args.gene:
        g = (r.get("genes") or {}).get(args.gene)
        if not g:
            print(f"{args.gene}: not a protein-coding gene in the dump (or unknown symbol)")
            return 1
        print(
            f"{args.gene}: origin {g['origin']} ({g['ladder']}); orthologues in {g['species']} species, "
            f"{g['one2one']} one-to-one; {len(g['paralogues'])} paralogues"
            + (f": {', '.join(g['paralogues'][: args.top])}" if g["paralogues"] else "")
            + (" …" if len(g["paralogues"]) > args.top else "")
        )
        print(f"  [{r['evidence']['origin']}; {r['evidence']['caveat']}]")
        return 0
    if args.library:
        lib = (r.get("libraries") or {}).get(args.library)
        if not lib:
            print(f"{args.library}: no members placed")
            return 1
        print(
            f"{args.library}: {lib['members_placed']} of {lib['members']} members placed; median origin "
            f"{lib['median_origin']} ({lib['ladder']}), deepest {lib['deepest']}; animal or older "
            f"{_pct(lib['share_animal_or_older'])}, vertebrate to tetrapod "
            f"{_pct(lib['share_vertebrate_to_tetrapod'])}, amniote or younger "
            f"{_pct(lib['share_amniote_or_younger'])}"
        )
        print(
            _table([{"stratum": s, "genes": n} for s, n in lib["distribution"].items()], ["stratum", "genes"])
        )
        return 0
    print(
        f"{r['coding_genes']:,} protein-coding genes placed by origin "
        f"({r['species_placed']} species on the ladder)"
    )
    print(
        _table(
            [
                {"stratum": s, "ladder": LADDER[s], "genes": r["by_stratum"].get(s, 0)}
                for s in STRATA
                if r["by_stratum"].get(s)
            ],
            ["stratum", "ladder", "genes"],
        )
    )
    rows = sorted(r["libraries"].items(), key=lambda kv: -kv[1]["share_animal_or_older"])
    print(f"\nlibraries by age (oldest first, top {args.top}):")
    print(
        _table(
            [
                {
                    "library": k,
                    "placed": v["members_placed"],
                    "median origin": v["median_origin"],
                    "animal+": _pct(v["share_animal_or_older"]),
                    "vertebrate": _pct(v["share_vertebrate_to_tetrapod"]),
                    "amniote+": _pct(v["share_amniote_or_younger"]),
                }
                for k, v in rows[: args.top]
            ],
            ["library", "placed", "median origin", "animal+", "vertebrate", "amniote+"],
        )
    )
    print(f"  [{r['evidence']['origin']}; {r['evidence']['caveat']}]")
    return 0


def cmd_variation(args: argparse.Namespace) -> int:
    """Constraint on two axes, mammals and people, per UNKNOWN block and element (area J step 1)."""
    from genomeos.attribution.variation import CASE_ORDER, run_and_save, vista_genome
    from genomeos.results import load_result, save_result

    if args.vista_genome:
        r = load_result("variation_vista_genome_wide") if args.cached else None
        if r is None:
            r = vista_genome(progress=lambda m: print(f"  {m}", flush=True))
            save_result("variation_vista_genome_wide", r)
        for k in ("positive", "negative"):
            a = r[k]
            print(
                f"VISTA {k}s: {a['intervals']} elements, {_pct(a['fraction_above'])} of their kilobases "
                f"constrained among people (Z >= {r['thresholds']['gnocchi']}), "
                f"{_pct(a['share_intervals_constrained'])} touch one, "
                f"{_pct(a['share_intervals_strong'])} a strong one"
            )
        print(f"  [{r['evidence']}; {r['cost']['mb_fetched']} MB in {r['cost']['seconds']} s]")
        return 0
    r = load_result(f"variation_{args.chrom}") if args.cached else None
    if r is None:
        r = run_and_save(args.chrom)
    print(
        f"{args.chrom}: {len(r['blocks'])} UNKNOWN blocks and {len(r['elements'])} attributed elements "
        f"read on both axes (mammals: phyloP >= {r['thresholds']['phylop']}; people: Gnocchi Z >= "
        f"{r['thresholds']['gnocchi']} over >= {_pct(r['thresholds']['human_min_fraction'], 0)} "
        "of the kilobases)"
    )
    rows = []
    for tier, v in r["by_tier"].items():
        cases = v.get("cases", {})
        rows.append(
            {
                "tier": tier,
                "blocks": v["blocks"],
                "Mb": f"{v['bp'] / 1e6:.2f}",
                "human-constrained": _pct(v["human_constrained_fraction"]),
                **{c: cases.get(c, {}).get("blocks", 0) for c in CASE_ORDER},
            }
        )
    print(_table(rows, ["tier", "blocks", "Mb", "human-constrained", *CASE_ORDER]))
    ctrl = r["controls"]
    if "canonical_cds" in ctrl:
        print(
            f"  controls: canonical coding segments {_pct(ctrl['canonical_cds']['fraction_above'])} of "
            f"kilobases constrained among people, introns "
            f"{_pct(ctrl['canonical_introns']['fraction_above'])}, "
            + ", ".join(f"{t} {_pct(f)}" for t, f in ctrl["unknown_by_tier"].items() if f is not None)
        )
    ec = r["by_element_case"]["by_case"]
    print(
        "  elements by case: "
        + "; ".join(
            f"{c} {ec[c]['elements']} (name a gene {_pct(ec[c]['name_a_gene_share'])})"
            for c in CASE_ORDER
            if ec[c]["elements"]
        )
        + f"; unmeasured {r['by_element_case']['unmeasured']}"
    )
    if r["vista"]:
        v = r["vista"]
        print(
            f"  VISTA on this chromosome: positives {_pct(v['positive']['fraction_above'])} of kilobases "
            f"constrained among people against negatives {_pct(v['negative']['fraction_above'])}"
        )
    for kv in r["known_values"]:
        print(f"  {kv['rsid']} ({kv['note']}): {kv['reading']}")
    top = sorted(
        (b for b in r["blocks"] if b["case"] and b["case"]["case"] in ("syntax", "relaxed")),
        key=lambda b: -(b["mammal_fraction"] or 0),
    )[: args.top]
    if top:
        print(f"\nmammal-constrained blocks by case (top {len(top)}):")
        print(
            _table(
                [
                    {
                        "locus": f"{args.chrom}:{b['start']:,}-{b['end']:,}",
                        "kb": f"{b['length'] / 1e3:.0f}",
                        "tier": b["tier"],
                        "mammals": _pct(b["mammal_fraction"]),
                        "people": _pct(b["gnocchi"]["fraction_above"]),
                        "max Z": b["gnocchi"]["maximum"],
                        "case": b["case"]["case"],
                        "conf": b["case"]["confidence"],
                    }
                    for b in top
                ],
                ["locus", "kb", "tier", "mammals", "people", "max Z", "case", "conf"],
            )
        )
    print(f"  [{r['evidence']['case']}; {r['cost']['seconds']} s]")
    return 0


def cmd_organise(args: argparse.Namespace) -> int:
    """The block organiser: every UNKNOWN block with both constraint axes and the copy flag, copies first."""
    from genomeos.attribution.organise import blocks, distil, run_and_save
    from genomeos.results import load_result, save_result

    if not args.chrom:
        s = distil()
        if not s["chromosomes"]:
            print("no organised_chr*.json yet: genomeos organise --chrom chr21")
            return 1
        save_result("organised_genome_wide", s)
        ru = s["real_unknown"]
        print(
            f"{s['chromosomes']} chromosomes organised; copies (>= {s['copy_min_fraction']:.0%} duplicated) "
            f"set apart; the real unknown after copies: {ru['blocks']:,} blocks, {ru['bp'] / 1e6:.1f} Mb"
        )
        rows = [
            {
                "tier": t,
                "blocks": v["blocks"],
                "copies": v["copies"],
                "after copies": v["blocks"] - v["copies"],
                **{
                    k: (v["cases_after_copies"].get(k) or {}).get("blocks", 0)
                    for k in ("syntax", "relaxed", "recent", "tolerant", "unmeasured")
                },
            }
            for t, v in s["by_tier"].items()
        ]
        print(
            _table(
                rows,
                [
                    "tier",
                    "blocks",
                    "copies",
                    "after copies",
                    "syntax",
                    "relaxed",
                    "recent",
                    "tolerant",
                    "unmeasured",
                ],
            )
        )
        return 0
    out = load_result(f"organised_{args.chrom}") if args.cached else None
    if out is None:
        out = run_and_save(args.chrom, top=args.top)
    ru = out["real_unknown"]
    print(
        f"{args.chrom}: {out['blocks']} UNKNOWN blocks, {out['with_human_axis']} read on the human axis, "
        f"{out['with_copy_flag']} with a copy flag; {out['copies']} copies ({out['copies_bp'] / 1e6:.2f} Mb) "
        f"set apart; the real unknown after copies: {ru['blocks']} blocks, {ru['bp'] / 1e3:.0f} kb"
    )
    rows = [
        {
            "tier": t,
            "blocks": v["blocks"],
            "copies": v["copies"],
            "after copies": v["after_copies"],
            **{
                k: (v["cases_after_copies"].get(k) or {}).get("blocks", 0)
                for k in ("syntax", "relaxed", "recent", "tolerant", "unmeasured")
            },
        }
        for t, v in out["by_tier"].items()
    ]
    print(
        _table(
            rows,
            [
                "tier",
                "blocks",
                "copies",
                "after copies",
                "syntax",
                "relaxed",
                "recent",
                "tolerant",
                "unmeasured",
            ],
        )
    )
    if out["candidates"]:
        print(f"\nthe real unknown, sharpest first (top {min(args.top, len(out['candidates']))}):")
        print(
            _table(
                [
                    {
                        "block": f"{r['start']:,}-{r['end']:,}",
                        "kb": f"{r['length'] / 1e3:.0f}",
                        "class": r["class"],
                        "mammals": _pct(r["mammal_fraction"]),
                        "people": _pct(r["human_fraction"]),
                        "case": r["case"] or "",
                        "reading": r["reading"][:70],
                    }
                    for r in out["candidates"][: args.top]
                ],
                ["block", "kb", "class", "mammals", "people", "case", "reading"],
            )
        )
    if args.blocks:
        print(f"\nevery block ({out['blocks']}):")
        print(
            _table(
                [
                    {
                        "block": f"{r['start']:,}-{r['end']:,}",
                        "tier": r["tier"],
                        "copy": "copy" if r["copy"] else "",
                        "case": r["case"] or "",
                        "elements": r["attributed_elements"] or "",
                        "reading": r["reading"][:60],
                    }
                    for r in blocks(args.chrom)
                ],
                ["block", "tier", "copy", "case", "elements", "reading"],
            )
        )
    print(
        "  [tier: the budget; case: Zoonomia against gnomAD Gnocchi (area J); copy: UCSC genomicSuperDups "
        "at >= 50% of the block; nothing recomputed; docs/ATTRIBUTION.md]"
    )
    return 0


def cmd_syntax(args: argparse.Namespace) -> int:
    """Syntax against values at one gene: constrained bases, where people differ, and the overlap."""
    from genomeos.attribution.syntax import save_path, syntax_values

    names = [n.strip() for n in args.name.split(",")] if args.name else None
    r = syntax_values(args.gene, args.chrom, names=names, flank=args.flank)
    if args.save:
        out = save_path(args.gene, r["people"])
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(r, indent=1))
        print(f"saved {out}")
    print(
        f"{r['gene']} {r['chrom']}:{r['span'][0]:,}-{r['span'][1]:,} ({r['bases'] / 1e3:.0f} kb, strand "
        f"{r['strand']}, flank {r['flank']}); people: {', '.join(r['people']) or 'none imported'}"
    )
    print(
        f"  syntax: {r['syntax_bases']:,} bases constrained across 241 mammals "
        f"({_pct(r['syntax_fraction'])}), mean phyloP {r['phylop_mean']}; "
        f"{r['features']['exons']} exons, {r['features']['ccres']} cCREs"
    )
    print(
        f"  values: {r['snv_positions']} SNV positions ({r['values_per_kb']} per kb) and "
        f"{r['indel_positions']} indels where these people differ from the reference; "
        f"{r['values_in_syntax']} of the SNVs sit on syntax ({_pct(r['syntax_hit_rate'])} against "
        f"{_pct(r['expected_if_random'])} if variation ignored syntax)"
    )
    rows = [
        {
            "pos": f"{x['pos']:,}",
            "ref>alt": f"{x['ref']}>{x['alt']}",
            "phyloP": x["phylop"] if x["phylop"] is not None else "",
            "class": x["class"],
            "where": x["where"],
            "carriers": ", ".join(f"{k} {v}" for k, v in x["carriers"].items()),
        }
        for x in r["rows"][: args.top]
    ]
    if rows:
        print(f"\nmost constrained variable positions (top {len(rows)}):")
        print(_table(rows, ["pos", "ref>alt", "phyloP", "class", "where", "carriers"]))
    for k in r["known_values"]:
        who = ", ".join(f"{a} {b}" for a, b in k["carriers"].items()) or "nobody here carries the alternative"
        where = (
            "inside the span" if k["inside_span"] else "outside this span (flank it, or ask the other gene)"
        )
        print(f"\n{k['id']} at {k['pos']:,}: {k['note']}; {where}; {who}")
    print(
        f"  [{r['evidence']['syntax']}; {r['evidence']['values']}; "
        f"{r['cost']['mb_fetched']} MB in {r['cost']['seconds']} s]"
    )
    return 0


def cmd_unknown(args: argparse.Namespace) -> int:
    from genomeos.genome import Annotation, Genome, default_gencode
    from genomeos.genome.unknown import investigate, load_patterns
    from genomeos.results import save_result

    genome = Genome.from_fasta(args.genome)
    chrom = args.chrom or next(iter(genome.chromosomes))
    gff = args.gff3 or default_gencode({chrom})
    ann = Annotation.from_gff3(gff, {chrom})
    patterns = load_patterns(args.patterns) if args.patterns else load_patterns()
    r = investigate(
        genome.chromosomes[chrom].sequence,
        ann,
        chrom,
        patterns,
        min_size=args.min_size,
        progress=lambda m: print(f"  {m}", flush=True),
        max_blocks=args.max_blocks,
    )
    print(f"{chrom}: {r['unknown_blocks']} UNKNOWN blocks, {r['unknown_bp']:,} bp")
    print(f"  classified {r['classified_fraction']:.1%}")
    print(
        _table(
            [
                {
                    "class": k,
                    "blocks": v["blocks"],
                    "bp": f"{v['bp']:,}",
                    "share": f"{v['bp'] / r['unknown_bp']:.1%}",
                }
                for k, v in r["by_class"].items()
            ],
            ["class", "blocks", "bp", "share"],
        )
    )
    print("largest blocks:")
    for b in r["blocks"][:10]:
        f = b["features"]
        where = f"{chrom}:{b['start']}-{b['end']}"
        feats = f"gc={f.get('gc')} tandem={f.get('tandem_fraction')} highcopy={f.get('high_copy_fraction')}"
        print(f"  {where:<24} {b['length']:>10,} bp  {b['class']:<22} conf={b['confidence']:.1f}")
        print(f"      {feats} orf={f.get('longest_orf_aa')}")
    name = args.save or f"unknown_{chrom}"
    print(f"  saved {save_result(name, r)}")
    return 0


def _print_definition(d: dict) -> None:
    s = d["sections"]
    ident = s["identity"]["items"] or {}
    print(
        f"{d['id'] or d['gene']}  {ident.get('name', '')}  {ident.get('length', '?')} aa  "
        f"[{s['identity']['evidence']}: {s['identity']['source']}]"
    )
    if ident.get("existence"):
        print(f"  existence: {ident['existence']}")
    go = (s.get("genomic_origin") or {}).get("items")
    if go:
        canon = next((t for t in go["transcripts"] if t["canonical"]), None)
        print(
            f"  origin: {go['gene_id']} {go['locus']}; {len(go['transcripts'])} transcripts, "
            f"{go['protein_products']} protein products; canonical {canon['name'] if canon else '?'} "
            f"[{s['genomic_origin']['evidence']}]"
        )
    iso = s.get("isoforms", {}).get("items") or []
    print(f"  isoforms (UniProt): {len(iso)}  {', '.join(i['id'] for i in iso[:6])}")
    fn = (s.get("function") or {}).get("items") or {}
    if fn.get("summary"):
        print(f"  function: {fn['summary'][0][:220]}…")
    if fn.get("location"):
        print(f"  location: {'; '.join(fn['location'][:4])}")
    dom = (s.get("domains") or {}).get("items") or {}
    print(
        f"  domains: {len(dom.get('interpro', []))} InterPro entries, {len(dom.get('features', []))} features"
        f"  [{s.get('domains', {}).get('evidence')}]"
    )
    print(f"  modifications: {len(s.get('modifications', {}).get('items') or [])} sites")
    exp = s.get("structures_experimental", {})
    pred = s.get("structures_predicted", {})
    methods: dict[str, int] = {}
    for x in exp.get("items") or []:
        methods[x["method"]] = methods.get(x["method"], 0) + 1
    print(
        f"  structures: experimental {exp.get('count', 0)} {dict(sorted(methods.items()))} "
        f"[experimental: PDB]; predicted {len(pred.get('items') or [])} [predicted: AlphaFold]"
    )
    pw = s.get("pathways", {}).get("items") or []
    print(f"  pathways (Reactome): {len(pw)}  {'; '.join(p['name'] for p in pw[:3])}")
    it = s.get("interactions", {})
    items = it.get("items") or []
    phys = sum(1 for x in items if x["physical_evidence"])
    print(
        f"  interactions (STRING ≥0.7): {len(items)}, {phys} with experimental support  "
        f"[{it.get('evidence')}]{'  (' + it['error'] + ')' if it.get('error') else ''}"
    )
    ex = s.get("expression", {})
    e = ex.get("items") or {}
    if e:
        print(
            f"  expression (HPA): {e.get('tissue_specificity')}; {e.get('tissue_distribution')}; "
            f"location {', '.join(e.get('subcellular_main') or []) or '?'}; "
            f"cell types: {e.get('cell_type_specificity')}"
        )
        if e.get("tissue_ntpm"):
            top = sorted(e["tissue_ntpm"].items(), key=lambda kv: -kv[1])[:5]
            print("    " + ", ".join(f"{k} {v:.0f}" for k, v in top))
    dis = s.get("diseases", {}).get("items") or []
    if dis:
        print(f"  diseases (UniProt): {len(dis)}  {'; '.join(x['name'] for x in dis[:4] if x['name'])}")
    cov = d.get("coverage", {})
    print("  coverage: " + "  ".join(f"{k} {'✓' if v else '·'}" for k, v in cov.items()))


def _accession_for(symbol_or_acc: str) -> str:
    import re

    if re.fullmatch(r"[A-NR-Z][0-9][A-Z0-9]{3}[0-9]([A-Z][A-Z0-9]{2}[0-9])?", symbol_or_acc):
        return symbol_or_acc
    from genomeos.molecules import compile_protein

    d = compile_protein(symbol_or_acc, sources={"uniprot"})
    ident = d["sections"].get("identity", {}).get("items")
    if not ident:
        raise SystemExit(f"no reviewed UniProt entry for {symbol_or_acc}")
    return ident["accession"]


def _pathway_kinetic(args: argparse.Namespace, ids: list[str]) -> int:
    """Kinetic counterpart: a curated BioModels ODE model for the pathway, run on the in-house engine."""
    from genomeos.lang.tools import spark
    from genomeos.molecules import biomodels
    from genomeos.molecules.reactome import PathwayModel, fetch_pathway
    from genomeos.results import save_result

    query = args.query
    if not query and ids:
        try:
            query = PathwayModel.from_sbml(fetch_pathway(ids[0])).summary()["name"]
            print(f"{ids[0]}: {query!r}")
        except Exception as e:  # noqa: BLE001
            print(f"{ids[0]}: unavailable ({str(e)[:80]}); give --query or --model")
            return 1
    model_id = args.model
    if not model_id:
        if not query:
            print("give a Reactome id, --query TEXT or --model BIOMD…")
            return 1
        hits, used = biomodels.search_pathway(query, limit=args.limit)
        if not hits:
            print(f"no curated BioModels entry matches {query!r}; try --query with other words")
            return 1
        print(f"BioModels, curated, matching {used!r}:")
        for h in hits:
            print(f"  {h['id']}  {h['name'][:80]}")
        model_id = hits[0]["id"]
        print(f"running the first ({model_id}); pick another with --model")
    path = biomodels.fetch(model_id)
    base = biomodels.run(path, duration=args.hours, dt=args.dt)
    print(
        f"{base['model']} {base['name']!r}: {base['species']} species, {base['reactions']} reactions, "
        f"{base['functions']} functions, {base['rate_rules']} rate rules; {args.hours} time units"
    )
    print(f"  [{base['evidence']}]")
    shown = sorted(base["levels"].items(), key=lambda kv: -kv[1]["peak"])[: args.top]
    for sid, lv in shown:
        name = lv["name"][:28]
        print(f"  {name:<28} {spark(base['series'][sid])}  final={lv['final']:10.3f}  peaks={lv['peaks']}")
    baseline = {k: v for k, v in base.items() if k not in ("series", "times")}
    out = {"pathways": ids, "query": query, "baseline": baseline}
    if args.knockout:
        acc = None
        try:
            acc = _accession_for(args.knockout)
        except (SystemExit, Exception):  # noqa: BLE001 - no compiled definition: match by name only
            acc = None
        ko = biomodels.run(
            path, duration=args.hours, dt=args.dt, knockout=[args.knockout], accessions={args.knockout: acc}
        )
        if not ko["knocked_out"]:
            print(f"  knockout {args.knockout}: no species of the model matches that name or accession {acc}")
        else:
            changed = biomodels.compare(base, ko)
            held = ", ".join(ko["knocked_out"])
            print(
                f"  knockout {args.knockout} → held at 0: {held}; {len(changed)} species change their final "
                f"level or peak by ≥ 10% of their range  [inferred from the run]"
            )
            for r in changed[: args.top]:
                name = r["name"][:28]
                print(
                    f"    {name:<28} final {r['baseline_final']:9.3f} → {r['knockout_final']:9.3f} "
                    f"({r['final_change']:+.0%})  peak {r['baseline_peak']:9.3f} → {r['knockout_peak']:9.3f} "
                    f"({r['peak_change']:+.0%})"
                )
            out["knockout"] = {"term": args.knockout, "held": ko["knocked_out"], "changed": changed}
    save_result(f"kinetic_{model_id}", out)
    return 0


def cmd_pathway(args: argparse.Namespace) -> int:
    from genomeos.molecules.reactome import PathwayModel, fetch_pathway
    from genomeos.results import save_result

    ids = list(args.ids)
    if args.gene:
        from genomeos.molecules import compile_protein

        d = compile_protein(args.gene)
        ids += [x["id"] for x in (d["sections"].get("pathways", {}).get("items") or [])][: args.limit]
        print(f"{args.gene}: {len(ids)} Reactome pathways from the compiled definition")
    if args.kinetic or args.model or args.query:
        return _pathway_kinetic(args, ids)
    if not ids:
        print("give pathway ids (R-HSA-…) or --gene SYMBOL")
        return 1
    acc = _accession_for(args.knockout) if args.knockout else None
    rows = []
    for pid in ids:
        try:
            model = PathwayModel.from_sbml(fetch_pathway(pid))
        except Exception as e:  # noqa: BLE001
            print(f"  {pid}: unavailable ({str(e)[:80]})")
            continue
        s = model.summary()
        if acc:
            k = model.knockout(acc)
            rows.append(
                {
                    "pathway": pid,
                    "name": s["name"][:50],
                    "reactions": s["reactions"],
                    "reachable": s["reactions_reachable_from_sources"],
                    "lost": len(k["reactions_lost"]),
                    "fraction_lost": f"{k['fraction_lost']:.0%}",
                    "unreachable_products": len(k["products_unreachable"]),
                }
            )
            if len(ids) == 1 or args.verbose:
                print(
                    f"{pid} {s['name']}: {s['species']} species ({s['proteins']} proteins), "
                    f"{s['reactions']} reactions, {s['reactions_reachable_from_sources']} reachable "
                    f"from sources  [{s['evidence']} v{s['version']}]"
                )
                print(
                    f"  knockout {acc}: {len(k['entities_containing'])} entities contain it; "
                    f"{len(k['reactions_lost'])} reactions lost ({k['fraction_lost']:.0%}), "
                    f"{len(k['products_unreachable'])} products unreachable  [{k['logic']}] "
                    f"conf={k['confidence']}"
                )
                for r in k["reactions_lost"][:15]:
                    print(f"    ✗ {r['name']}  ({r['id']})")
                for sp in k["products_unreachable"][:10]:
                    print(f"    · {sp}")
        else:
            rows.append(
                {
                    "pathway": pid,
                    "name": s["name"][:50],
                    "species": s["species"],
                    "proteins": s["proteins"],
                    "reactions": s["reactions"],
                    "reachable": s["reactions_reachable_from_sources"],
                    "catalysed": s["catalysed"],
                    "inhibited": s["inhibited"],
                }
            )
    if len(rows) > 1:
        key = "lost" if acc else "reactions"
        rows.sort(key=lambda r: -r[key])
        print(_table(rows, list(rows[0].keys())))
    if args.gene and acc:
        save_result(
            f"knockout_{args.gene.upper()}",
            {
                "gene": args.gene.upper(),
                "accession": acc,
                "pathways": rows,
                "evidence": "curated: Reactome; logic inferred (reachability)",
            },
        )
        print(f"  saved data/results/knockout_{args.gene.upper()}.json")
    return 0


def cmd_proteome(args: argparse.Namespace) -> int:
    from genomeos.molecules.proteome import QUESTIONS, compile_chromosome
    from genomeos.results import save_result

    r = compile_chromosome(args.chrom, args.gff3, args.limit)
    print(
        f"{args.chrom}: {r['coding_genes']} coding genes compiled in {r['seconds']} s; "
        f"{len(r['no_reviewed_entry'])} without a reviewed UniProt entry"
    )
    print(
        _table(
            [
                {
                    "question": q,
                    "proteins": r["coverage_counts"][q],
                    "fraction": f"{r['coverage_fraction'][q]:.1%}"
                    if r["coverage_fraction"][q] is not None
                    else "-",
                }
                for q in QUESTIONS
            ],
            ["question", "proteins", "fraction"],
        )
    )
    if not args.limit:
        save_result(f"proteome_{args.chrom}", r)
        print(f"  saved data/results/proteome_{args.chrom}.json")
    return 0


def cmd_protein(args: argparse.Namespace) -> int:
    from genomeos.genome import Annotation, IndexedGenome, default_gencode
    from genomeos.molecules import protein_report

    if args.lib:
        from genomeos.lib.proteome import block, get

        r = get(args.symbol)
        if not r:
            print(f"{args.symbol.upper()}: not in the packaged proteome (19,283 human proteins)")
            return 1
        if args.bio:
            print(block(args.symbol), end="")
            return 0
        print(
            f"{args.symbol.upper()}  UniProt:{r['accession']}  {r['name']}  {r['length']} aa  [{r['chrom']}]"
        )
        print(f"  {r['function']}")
        print(
            f"  location: {'; '.join(r['location']) or '?'}  · tissue: {r['tissue_pattern'] or '?'}  · "
            f"cell types: {r['cell_type_pattern'] or '?'}"
        )
        print(f"  domains: {', '.join(r['domains']) or 'none'}")
        print(
            f"  pathways: {len(r['pathways'])}  · partners (physical): {', '.join(r['partners']) or 'none'}"
        )
        af = ", AlphaFold" if r["alphafold"] else ""
        print(
            f"  structures: {r['structures_experimental']} experimental{af}  · "
            f"diseases: {'; '.join(r['diseases']) or 'none recorded'}"
        )
        print("  [packaged proteome: UniProt, InterPro, Reactome, STRING, HPA, PDB, AlphaFold; no network]")
        return 0
    if args.bio:
        from genomeos.molecules import compile_protein, to_biolang

        print(to_biolang(compile_protein(args.symbol, refresh=args.refresh)), end="")
        return 0
    if args.compile:
        from genomeos.molecules import compile_protein, states_from_definition

        d = compile_protein(args.symbol, refresh=args.refresh)
        _print_definition(d)
        st = states_from_definition(d)
        if st:
            print(
                f"  {len(st)} ProteinState records derivable (tissue / cell type levels), "
                f"e.g. {st[0].to_dict()}"
            )
        return 0

    ann = genome = None
    if args.chrom:
        gff = default_gencode({args.chrom})
        if gff:
            ann = Annotation.from_gff3(gff, {args.chrom})
            genome = IndexedGenome(args.genome)
    r = protein_report(args.symbol, ann, genome, structure=not args.no_structure)
    if genome:
        genome.close()
    if r["ours"]:
        o = r["ours"]
        print(f"{args.symbol}: our translation of {o['transcript']}: {o['length']} aa")
        print(f"  ({o['coding_transcripts']} coding transcripts annotated)")
    u = r["uniprot"]
    if not u:
        print("  no reviewed UniProt entry")
        return 1
    print(f"  UniProt {u['accession']} {u['name']}: {u['length']} aa  [{u['evidence']}]")
    if r.get("agreement"):
        a = r["agreement"]
        print(
            f"  agreement with our translation: identity {a['identity']:.1%}, same length {a['same_length']}"
        )
    if u["location"]:
        print(f"  location: {u['location'][:120]}")
    if u["function"]:
        print(f"  function: {u['function'][:240]}…")
    for f in u["features"][:12]:
        print(f"    {f['type']:<20} {f['start']:>5}-{f['end']:<5} {f['description'][:50]}")
    s = r.get("structure")
    if s:
        print(
            f"  AlphaFold {s['entry']} ({s['version']}): {s['length']} residues, mean pLDDT {s['mean_plddt']}"
        )
        print(f"    {s['confident_fraction']:.0%} confident (>=70), Rg {s['radius_of_gyration_A']} Å")
        print("    [predicted]")
    return 0


def cmd_flow(args: argparse.Namespace) -> int:
    """Walk one gene upward: sequence → regulation → RNA → protein → cell → tissue → organism."""
    from genomeos.coords import Locus
    from genomeos.genome import Annotation, IndexedGenome, default_gencode
    from genomeos.genome.regulatory import ccre_index, count_in, load_ccres
    from genomeos.lib import LIBRARIES, KnowledgeBase
    from genomeos.runtime import coding_sequence, splice, translate

    gff = args.gff3 or default_gencode({args.chrom})
    if not gff:
        print("no annotation for that chromosome locally")
        return 1
    ann = Annotation.from_gff3(gff, {args.chrom})
    genome = IndexedGenome(args.genome)
    try:
        g = ann.gene(args.symbol)
    except KeyError:
        print(f"{args.symbol} not on {args.chrom}")
        return 1
    m = ann.to_module("flow")
    txs = [t for t in m.entities[g.id].transcripts if t.cds_segments]
    canon = next((t for t in txs if "Ensembl_canonical" in t.tags), txs[0] if txs else None)
    print(f"1 DNA        {g.symbol} {g.locus}  {g.locus.length:,} bp  [curated: GENCODE]")
    head = str(genome.fetch(Locus(args.chrom, g.locus.start, g.locus.start + 60)))
    print(f"             {head}…")
    # regulation: ENCODE elements around the TSS and inside the gene; CpG island at the promoter
    tss = g.locus.end if g.locus.strand.value == "-" else g.locus.start
    idx = ccre_index(load_ccres(args.chrom))
    near = count_in(idx, tss - 2000, tss + 2000) if idx else {}
    inside = count_in(idx, g.locus.start, g.locus.end) if idx else {}
    prom = str(genome.fetch(Locus(args.chrom, max(0, tss - 500), tss + 500)))
    c_, g_, cg = prom.count("C"), prom.count("G"), prom.count("CG")
    island = c_ and g_ and (c_ + g_) / len(prom) > 0.5 and cg * len(prom) / (c_ * g_) > 0.6
    print(f"2 regulation promoter CpG island: {'yes' if island else 'no'} [predicted]")
    print(
        f"             ENCODE elements ±2 kb of TSS: {near or 'none'}; inside the gene: "
        f"{sum(inside.values()) if inside else 0} [curated: ENCODE]"
    )
    if idx:
        from genomeos.genome.regulation import regulation_of

        reg = regulation_of(g.symbol, args.chrom, load_ccres(args.chrom), ann, genome.lengths[args.chrom])
        dom = reg["domain"]
        print(
            f"             node {dom['id'] if dom else '?'}: {reg['enhancers_in_domain']} enhancers can "
            f"reach it ({reg['enhancers_nearest_to_this_gene']} nearest to this gene), "
            f"{len(reg['promoters'])} promoter elements, {len(reg['insulators_bounding'])} bounding "
            f"insulators [inferred: reach bounded by CTCF domain]"
        )
    # RNA
    if canon is None:
        print("3 RNA        no coding transcript")
        genome.close()
        return 0
    mrna = splice(genome, canon)
    cds_nt = len(coding_sequence(genome, canon))
    print(
        f"3 RNA        {len(g.transcripts)} transcripts; canonical {canon.attrs.get('name', canon.id)}: "
        f"{len(canon.exons)} exons"
    )
    print(f"             spliced to {len(mrna):,} nt mRNA, CDS {cds_nt:,} nt  [curated: GENCODE]")
    # protein
    prot = translate(coding_sequence(genome, canon), initiator=True)
    print(f"4 protein    {len(prot)} aa: {prot[:40]}…  [translated by GenomeOS from the DNA above]")
    try:
        from genomeos.molecules import uniprot_entry

        u = uniprot_entry(g.symbol)
        if u:
            ident = sum(1 for a, b in zip(prot, u.sequence, strict=False) if a == b) / max(
                len(prot), len(u.length and u.sequence)
            )
            print(
                f"             UniProt {u.accession} {u.name}: {u.length} aa, identity {ident:.1%}  [curated]"
            )
    except Exception as e:  # noqa: BLE001
        print(f"             (UniProt unavailable: {e})")
    # cell behaviour: libraries
    kb = KnowledgeBase() if KnowledgeBase.available() else None
    libs = kb.libraries_of(g.symbol) if kb else []
    layers = {}
    for lib_id in libs:
        layers.setdefault(LIBRARIES[lib_id].layer, []).append(lib_id.split(".", 1)[1])
    print(f"5 cell       libraries (GO/Reactome by data): core {', '.join(layers.get('core', [])) or '-'}")
    print(f"             timer {', '.join(layers.get('timer', [])) or '-'}")
    print(f"6 tissue     blueprint {', '.join(layers.get('blueprint', [])) or '-'}")
    print(f"             systems {', '.join(layers.get('systems', [])) or '-'}")
    print("7 organism   see `genomeos twin` (a person's state), `genomeos organism` (lineage),")
    print("             `genomeos anatomy --compare` (the whole genome's organisation)")
    genome.close()
    return 0


def cmd_segments(args: argparse.Namespace) -> int:
    """Segment parser v1: chain learned signals into candidate genes and score them against GENCODE."""
    from genomeos.genome import Annotation, IndexedGenome, default_gencode
    from genomeos.genome.segments import parse_chromosome
    from genomeos.genome.signals import SignalSet
    from genomeos.results import save_result

    gff = default_gencode({args.chrom})
    if not gff or not Path(args.genome).exists() or not Path(args.signals).exists():
        print(f"{args.chrom}: needs local models, sequence and learned signals ({args.signals})")
        return 1
    ann = Annotation.from_gff3(gff, {args.chrom})
    g = IndexedGenome(args.genome)
    from genomeos.coords import Locus

    length = g.lengths[args.chrom]
    seq = str(g.fetch(Locus(args.chrom, 0, length)))
    g.close()
    sig = SignalSet.load(args.signals)
    sites = None
    name = f"segments_{args.chrom}"
    if args.predicted_sites:
        from genomeos.predict import status
        from genomeos.predict.splice_sites import SpliceSites, client_factory

        st = status()
        if not st["enabled"]:
            print(f"AlphaGenome splice sites are disabled: {st['reason']}. To enable: {st['how']}")
            return 2
        ss = SpliceSites(args.chrom, length).load(
            client_factory, progress=lambda m: print("  " + m, flush=True)
        )
        sm = ss.summary()
        print(
            f"{args.chrom}: predicted splice sites over {sm['windows']} windows "
            f"({sm['requests_this_run']} requests, {sm['seconds_this_run']} s); "
            f"confident (p ≥ 0.5): {sm['confident_at_or_above_0.5']}"
        )
        sites = ss.local
        name = f"segments_{args.chrom}_predicted_sites"
    start_windows = None
    if args.promoter_anchored is not None:
        from genomeos.genome.regulatory import load_ccres
        from genomeos.genome.segments import promoter_start_windows

        ccres = load_ccres(args.chrom)
        if not ccres:
            print(f"{args.chrom}: no ENCODE elements (genomeos data fetch --chrom {args.chrom})")
            return 1
        start_windows = promoter_start_windows(ccres, args.promoter_anchored)
        allowed = sum(b - a for a, b in start_windows["+"])
        print(
            f"{args.chrom}: gene starts confined to {allowed / length:.1%} of the chromosome "
            f"(ENCODE promoter-like elements, {args.promoter_anchored:,} bp downstream)"
        )
        name += "_anchored"
    post_filter = None
    if args.rna_filter is not None:
        from genomeos.predict import status
        from genomeos.predict.rna_tracks import RnaCoverage, filter_predictions
        from genomeos.predict.splice_sites import client_factory

        st = status()
        if not st["enabled"]:
            print(f"AlphaGenome RNA tracks are disabled: {st['reason']}. To enable: {st['how']}")
            return 2
        cov = RnaCoverage(args.chrom, length).load(
            client_factory, progress=lambda m: print("  " + m, flush=True)
        )
        cs = cov.summary()
        print(
            f"{args.chrom}: predicted RNA over {cs['windows']} windows ({cs['requests_this_run']} requests); "
            f"transcribed bases + {cs['transcribed_bases']['+']:,}, − {cs['transcribed_bases']['-']:,} "
            f"(coverage ≥ {cs['coverage_threshold']} in a panel of {len(cs['panel'])} tissues)"
        )
        post_filter = lambda preds: filter_predictions(preds, cov, args.rna_filter)  # noqa: E731
    measured = None
    if args.rna_measured:
        from genomeos.genome.rna_measured import MeasuredRna

        measured = MeasuredRna(args.rna_measured, args.chrom)
        frac = args.rna_filter if args.rna_filter is not None else 0.3
        inner = post_filter

        def post_filter(preds, measured=measured, frac=frac, inner=inner):  # noqa: F811
            if inner is None:
                return measured.filter(preds, frac, progress=lambda d, n, b: None)
            predicted_keep = inner(preds)
            if args.rna_combine == "intersection":
                return measured.filter(predicted_keep, frac, progress=lambda d, n, b: None)
            # union: the model's gene-level tracks for precision, the measured lines for the genes the panel
            # does not express; a candidate stays if either says it is transcribed
            keep_ids = {id(p) for p in predicted_keep}
            measured_keep = measured.filter(preds, frac, progress=lambda d, n, b: None)
            keep_ids |= {id(p) for p in measured_keep}
            return [p for p in preds if id(p) in keep_ids]

    r = parse_chromosome(
        args.chrom,
        seq,
        sig,
        ann,
        min_relative=args.min_relative,
        sites=sites,
        start_windows=start_windows,
        coding_model=args.coding,
        post_filter=post_filter,
    )
    if measured is not None:
        ms = measured.summary()
        print(
            f"{args.chrom}: measured RNA of {ms['cell_type']} over {ms['exons_measured']:,} candidate exons "
            f"({ms['bytes_fetched'] / 1e6:.1f} MB of bigWig read; tracks {ms['tracks']})"
        )
    if args.coding != "codon":
        name += f"_{args.coding}"
    if args.rna_filter is not None:
        name += "_rna"
        r["rna"] = {**cs, "min_exon_fraction": args.rna_filter}
        r["evidence"] += "; candidates kept only where AlphaGenome predicts RNA over their exons (predicted)"
    if measured is not None:
        cells = [x.strip() for x in args.rna_measured.split(",") if x.strip()]
        name += "_measured_" + (cells[0].replace(" ", "_") if len(cells) == 1 else f"panel{len(cells)}")
        if args.rna_filter is not None:
            name += f"_{args.rna_combine}"
            r["rna_combine"] = args.rna_combine
        r["rna_measured"] = {**ms, "min_exon_fraction": frac}
        r["evidence"] += f"; candidates kept only where ENCODE {args.rna_measured} RNA-seq covers their exons"
    if args.predicted_sites:
        r["splice_sites"] = sm
    if start_windows is not None:
        r["promoter_anchored_bp"] = args.promoter_anchored
        r["evidence"] += "; gene starts confined to ENCODE promoter-like elements (curated)"
    save_result(name, r)
    print(
        f"{args.chrom}: {r['predictions']} candidate genes from signals alone, {r['predicted_exons']} exons"
    )
    print(
        f"  exons: sensitivity {r['exon_sensitivity']:.1%} precision {r['exon_precision']:.1%} "
        f"(of {r['true_exons']} annotated CDS segments)"
    )
    print(f"  splice sites: sensitivity {r['site_sensitivity']:.1%} precision {r['site_precision']:.1%}")
    print(
        f"  genes: sensitivity {r['gene_sensitivity']:.1%} precision {r['gene_precision']:.1%} "
        f"(of {r['true_genes']} coding genes)"
    )
    print(f"  [{r['evidence']}]  saved data/results/{name}.json")
    return 0


def cmd_reader(args: argparse.Namespace) -> int:
    """Reader v1: which nodes and genes a cell type reads, from ENCODE DNase peaks."""
    from genomeos.genome import Annotation, IndexedGenome, default_gencode
    from genomeos.genome.domains import infer_domains
    from genomeos.genome.reader import compare, fetch_peaks, load_peaks, read_chromosome, slug
    from genomeos.genome.regulatory import load_ccres
    from genomeos.results import save_result

    gff = default_gencode({args.chrom})
    if not gff or not Path(args.genome).exists():
        print(f"{args.chrom}: needs local models and sequence (genomeos data fetch --chrom {args.chrom})")
        return 1
    ann = Annotation.from_gff3(gff, {args.chrom})
    g = IndexedGenome(args.genome)
    length = g.lengths[args.chrom]
    g.close()
    ccres = load_ccres(args.chrom)
    domains = infer_domains(args.chrom, length, ccres, ann) if ccres else []
    results = []
    for cell in [args.cell_type, *(args.versus or [])]:
        if not load_peaks(cell, args.chrom):
            print(f"fetching ENCODE DNase peaks for {cell}…", flush=True)
            info = fetch_peaks(cell, {args.chrom})
            kept = info["kept"][args.chrom]
            print(f"  {info['accession']} ({info['candidates']} candidate files); kept {kept:,} peaks")
        r = read_chromosome(cell, args.chrom, ann, domains, ccres)
        save_result(f"reader_{slug(cell)}_{args.chrom}", {k: v for k, v in r.items() if k != "_read_all"})
        genes = f"{r['genes_read']}/{r['coding_genes']} coding genes read ({r['read_fraction']:.0%})"
        frac = f"{r['enhancers_active_fraction']:.0%}"
        enh = f"enhancers active {r['enhancers_active']}/{r['enhancers']} ({frac})"
        nodes = f"nodes open {r['nodes_open']}/{r['nodes']}, silent {r['nodes_silent']}"
        print(f"{cell} on {args.chrom}: {r['peaks']:,} peaks; {genes}; {enh}; {nodes}")
        print(
            "  strongest promoters: "
            + ", ".join(f"{x['gene']} {x['signal']:.0f}" for x in r["top_read"][:10])
        )
        print(f"  [{r['evidence']['peaks']}; {r['evidence']['read']}]")
        results.append(r)
    for other in results[1:]:
        c = compare(results[0], other)
        print(
            f"{c['a']} vs {c['b']}: {c['read_in_both']} genes read in both; "
            f"read only in {c['a']}: {', '.join(c['read_in_a_only'][:25])}"
        )
        print(f"  read only in {c['b']}: {', '.join(c['read_in_b_only'][:25])}")
        save_result(f"reader_{slug(c['a'])}_vs_{slug(c['b'])}_{args.chrom}", c)
    return 0


def cmd_report(args: argparse.Namespace) -> int:
    """Every layer about one gene, as one Markdown dossier."""
    from genomeos.genome import Annotation, IndexedGenome, default_gencode
    from genomeos.report import gene_report, to_markdown

    gff = default_gencode({args.chrom})
    if not gff or not Path(args.genome).exists():
        print(f"{args.chrom}: needs local models and sequence (genomeos data fetch --chrom {args.chrom})")
        return 1
    ann = Annotation.from_gff3(gff, {args.chrom})
    genome = IndexedGenome(args.genome)
    try:
        rep = gene_report(args.symbol, args.chrom, ann, genome)
    except KeyError:
        print(f"{args.symbol} not on {args.chrom}")
        return 1
    finally:
        genome.close()
    text = json.dumps(rep, indent=1) if args.json else to_markdown(rep)
    if args.out:
        Path(args.out).write_text(text)
        print(f"wrote {args.out}")
    else:
        print(text)
    return 0


def cmd_verify(args: argparse.Namespace) -> int:
    """Translate every compiled gene of a chromosome and compare with UniProt."""
    from genomeos.genome import Annotation, IndexedGenome, default_gencode
    from genomeos.molecules.verify import verify_chromosome
    from genomeos.results import save_result

    gff = default_gencode({args.chrom})
    if not gff or not Path(args.genome).exists():
        print(f"{args.chrom}: needs local models and sequence (genomeos data fetch --chrom {args.chrom})")
        return 1
    ann = Annotation.from_gff3(gff, {args.chrom})
    genome = IndexedGenome(args.genome)
    r = verify_chromosome(args.chrom, ann, genome)
    genome.close()
    save_result(f"translation_vs_uniprot_{args.chrom}", r)
    print(
        f"{args.chrom}: {r['genes_checked']} compiled genes; canonical transcript identical to UniProt "
        f"{r['exact_canonical_fraction']:.1%}; some isoform identical "
        f"{r['exact_some_isoform_fraction']:.1%}; "
        f"{r['canonical_differs_but_another_isoform_matches']} canonical-choice differences; "
        f"{r['same_protein_different_boundaries']} same protein with other boundaries; "
        f"{r['disagreement_count']} real disagreements  [{r['evidence']}]"
    )
    if r["disagreements"]:
        print(
            _table(
                [
                    {
                        "gene": d["gene"],
                        "transcript": d["transcript"],
                        "ours aa": d["ours_aa"],
                        "UniProt aa": d["uniprot_aa"],
                        "best identity": f"{d['identity_best_isoform']:.1%}",
                        "similarity": f"{d['similarity_best_isoform']:.1%}",
                    }
                    for d in r["disagreements"][: args.top]
                ],
                ["gene", "transcript", "ours aa", "UniProt aa", "best identity", "similarity"],
            )
        )
    print(f"  saved data/results/translation_vs_uniprot_{args.chrom}.json")
    return 0


def cmd_individual(args: argparse.Namespace) -> int:
    """A person's genome as a local individual: import once, then every layer sees it (offline)."""
    from genomeos.genome import individuals as ind

    if args.action == "import":
        try:
            m = ind.import_vcf(
                args.vcf,
                args.name,
                args.note or "",
                replace=args.replace,
                progress=lambda s: print(s, flush=True),
            )
        except (ValueError, FileNotFoundError, FileExistsError) as ex:
            print(ex)
            return 1
        chroms = ", ".join(f"{k} {v:,}" for k, v in m["chromosomes"].items())
        print(f"imported {m['name']}: {m['variants']:,} PASS variants on {len(m['chromosomes'])} chromosomes")
        print(f"  {chroms}")
        print(
            f"  skipped {m['skipped_filtered']:,} filtered rows and {m['skipped_other_contigs']:,} rows "
            f"on other contigs; {m['phased']:,} phased genotypes; sample column {m['sample_column']}"
        )
        print(f"  kept under data/individuals/{m['name']} (git-ignored; nothing leaves this machine)")
        print("  now: genomeos lookup, report, twin build --vcf, and individual genes see this person")
        return 0
    if args.action == "remove":
        ok = ind.remove(args.name)
        print(f"removed {args.name}" if ok else f"{args.name}: not an imported individual")
        return 0 if ok else 1
    if args.action == "predict":
        return _individual_predict(args)
    if args.action == "report":
        try:
            md = ind.dossier(args.name)
        except FileNotFoundError as ex:
            print(ex)
            return 1
        if args.out:
            Path(args.out).write_text(md)
            print(f"wrote {args.out}")
        else:
            print(md)
        return 0
    if args.action == "diff":
        from genomeos.genome.diff import genome_diff, render

        try:
            d = genome_diff(args.name, args.against, progress=lambda m: print("  " + m, flush=True))
        except FileNotFoundError as ex:
            print(ex)
            return 1
        md = render(d, top=args.top)
        if args.out:
            Path(args.out).write_text(md)
            print(f"wrote {args.out}")
        else:
            print(md)
        return 0
    if args.action == "protein":
        try:
            r = ind.tissue_proteins(args.name, args.gene, args.chrom)
        except (FileNotFoundError, KeyError) as ex:
            print(ex)
            return 1
        print(
            f"{r['individual']} × {r['gene']}: the protein each of {r['tissues_measured']} tissues makes  "
            f"[{r['evidence']}]"
        )
        for p in r["proteins"]:
            if "note" in p:
                print(f"  {p['transcript']}: {p['note']} (dominant in {len(p['tissues'])} tissues)")
                continue
            tag = " canonical" if p["canonical"] else ""
            ch = (
                "; ".join(
                    f"{v['hgvs_p'] or v['consequence']} ({v['genotype']})" for v in p["protein_changing"]
                )
                or "none"
            )
            print(
                f"  {p['name'] or p['transcript']}{tag}: {p['protein_length']} aa, dominant in "
                f"{p['tissues_dominant']} tissues ({', '.join(p['tissues'][:4])}"
                f"{'…' if len(p['tissues']) > 4 else ''}); "
                f"{p['coding_snvs']} coding SNVs, protein-changing: {ch}"
            )
        print(f"  {r['note']}")
        return 0
    if args.action == "regions":
        try:
            counts = ind.import_regions(args.name, args.bed, progress=lambda m: print(m, flush=True))
        except FileNotFoundError as ex:
            print(ex)
            return 1
        print(
            f"  {sum(counts.values()):,} intervals kept under "
            f"data/individuals/{args.name}/regions_<chrom>.bed"
        )
        return 0
    if args.action == "trio":
        try:
            r = ind.trio(args.name, args.father, args.mother, args.chrom or None)
        except FileNotFoundError as ex:
            print(ex)
            return 1
        t = r["totals"]
        print(
            f"{r['child']} with {r['father']} and {r['mother']} on {len(r['chromosomes'])} chromosomes: "
            f"{t['child_variants']:,} child variants; {r['fraction_inherited']:.1%} inherited "
            f"({t['in_both_parents']:,} from both), {t['de_novo_candidates']:,} absent from both parents "
            f"({r['fraction_de_novo_candidates']:.2%}), {t['mendelian_errors']:,} Mendelian errors "
            f"({r['fraction_mendelian_errors']:.3%}); {t['outside_a_parent_region']:,} outside a parent's "
            f"trusted region  [{r['regions']}]"
        )
        if t.get("de_novo_snv") is not None:
            print(
                f"  de novo candidates: {t['de_novo_snv']:,} SNVs, {t['de_novo_indel']:,} indels  "
                f"[{r.get('representation', '')}]"
            )
        rows = [
            {
                "chr": c,
                "child": f"{v['child_variants']:,}",
                "inherited": f"{v['inherited']:,}",
                "both": f"{v['in_both_parents']:,}",
                "de novo?": v["de_novo_candidates"],
                "errors": v["mendelian_errors"],
                "untrusted": v["outside_a_parent_region"],
            }
            for c, v in list(r["per_chromosome"].items())[: args.top]
        ]
        print(_table(rows, ["chr", "child", "inherited", "both", "de novo?", "errors", "untrusted"]))
        print(f"  [{r['evidence']}]  {r['note']}")
        print(f"  stored under data/individuals/{r['child']}/trio_{r['father']}_{r['mother']}.json")
        return 0
    if args.action == "coding":
        try:
            r = ind.coding_inventory(
                args.name,
                args.chrom or None,
                predict=getattr(args, "predict", False),
                progress=lambda m: print("  " + m, flush=True),
            )
        except FileNotFoundError as ex:
            print(ex)
            return 1
        bc = r["by_consequence"]
        print(
            f"{r['individual']}: {sum(bc.values()):,} coding SNVs on canonical transcripts over "
            f"{len(r['chromosomes'])} chromosomes; {r['genes_with_protein_changing']:,} genes with a "
            f"protein-changing variant, {r['genes_with_homozygous_changing']:,} with a homozygous one"
        )
        print("  by consequence: " + ", ".join(f"{k} {v:,}" for k, v in bc.items()))
        rows = [
            {
                "gene": x["gene"],
                "chrom": x["chrom"],
                "changing": x["protein_changing"],
                "homozygous": x["homozygous_changing"],
                "coding": x["coding_snvs"],
                "aa": x["protein_length"],
                "examples": "; ".join(x["examples"][:3])[:60],
            }
            for x in r["top"][: args.top]
        ]
        print(_table(rows, ["gene", "chrom", "changing", "homozygous", "coding", "aa", "examples"]))
        if r.get("missense_by_effect"):
            c = r["missense_predicted"]
            print(
                f"  missense predicted (AlphaMissense): {c['likely_pathogenic']:,} likely pathogenic "
                f"({r['missense_likely_pathogenic_homozygous']:,} homozygous), {c['ambiguous']:,} ambiguous, "
                f"{c['likely_benign']:,} likely benign, {c['unscored']:,} unscored; strongest first:"
            )
            print(
                _table(
                    [
                        {
                            "gene": m["gene"],
                            "change": m["hgvs_p"] or (m["predicted"] or {}).get("protein_variant") or "",
                            "genotype": m["genotype"],
                            "score": f"{m['predicted']['score']:.3f}" if m.get("predicted") else "",
                            "class": ((m.get("predicted") or {}).get("class") or "unscored").replace(
                                "_", " "
                            ),
                            "site": m["site"],
                        }
                        for m in r["missense_by_effect"][: args.top]
                    ],
                    ["gene", "change", "genotype", "score", "class", "site"],
                )
            )
        mr = [m for m in r.get("missense_ranked", []) if m["site"] != "none"][: args.top]
        if mr:
            print(
                f"  missense on an annotated UniProt site: {r['missense_at_annotated_site']:,} of "
                f"{r['missense']:,} ({r['missense_in_domain']:,} more inside a domain); read these first:"
            )
            print(
                _table(
                    [
                        {
                            "gene": m["gene"],
                            "change": m["hgvs_p"] or "",
                            "gt": m["genotype"],
                            "where": m["site"],
                            "annotation": "; ".join(m["features"][:2])[:70],
                        }
                        for m in mr
                    ],
                    ["gene", "change", "gt", "where", "annotation"],
                )
            )
        rt = r.get("reference_truncating")
        if rt:
            hidden = [x["gene"] for x in r.get("reference_truncating_rows", []) if x["matches_reference"]]
            print(
                f"  reference truncating alleles: hg38 carries one in {rt['genes']} genes; this person "
                f"matches the reference in {rt['matches_reference']} of {rt['checked']} on file "
                f"(carries the truncation): " + ", ".join(hidden[:14]) + ("…" if len(hidden) > 14 else "")
            )
        print(f"  [{r['evidence']}]  {r['note']}")
        return 0
    if args.action == "knockouts":
        try:
            r = ind.knockouts(args.name, args.chrom or None)
        except FileNotFoundError as ex:
            print(ex)
            return 1
        print(
            f"{r['individual']}: {len(r['hits'])} truncating SNVs on canonical transcripts over "
            f"{len(r['chromosomes'])} chromosomes ({r['nonsense']} nonsense, {r['start_lost']} start lost, "
            f"{r['stop_lost']} stop lost; {r['homozygous']} homozygous) in {len(r['genes'])} genes"
        )
        rows = [
            {
                "gene": h["gene"],
                "variant": f"{h['chrom']}:{h['pos']:,} {h['ref']}>{h['alt']}",
                "gt": h["genotype"],
                "zygosity": h["zygosity"],
                "consequence": h["consequence"],
                "protein": f"{h['hgvs_p'] or ''} of {h['protein_length']} aa",
                "lost": f"{h['fraction_lost']:.0%}" if h["fraction_lost"] is not None else "",
            }
            for h in r["hits"][: args.top]
        ]
        print(_table(rows, ["gene", "variant", "gt", "zygosity", "consequence", "protein", "lost"]))
        print(f"  [{r['evidence']}]  {r['note']}")
        print(f"  stored under data/individuals/{args.name}/knockouts.json (never under data/results)")
        return 0
    if args.action == "screen":
        from genomeos.genome import clinvar

        if not clinvar.pathogenic_path().exists():
            print("ClinVar not distilled yet: streaming the GRCh38 VCF once (about 190 MB, kept locally)")
            s = clinvar.distil(progress=lambda m: print("  " + m, flush=True))
            print(f"  kept {s['rows_kept']:,} pathogenic / likely pathogenic rows of {s['rows_read']:,}")
        try:
            r = clinvar.screen(args.name, args.chrom or None)
        except FileNotFoundError as ex:
            print(ex)
            return 1
        print(
            f"{r['individual']}: {len(r['hits'])} ClinVar pathogenic alleles carried over "
            f"{r['variants_scanned']:,} variants on {len(r['chromosomes'])} chromosomes "
            f"({r['pathogenic']} pathogenic, {r['likely_pathogenic']} likely; {r['homozygous']} homozygous; "
            f"{r['two_stars_or_more']} with ≥ 2 stars)"
        )
        rows = [
            {
                "variant": f"{h['chrom']}:{h['pos']:,} {h['ref'][:8]}>{h['alt'][:8]}",
                "gene": h["gene"],
                "gt": h["genotype"],
                "zygosity": h["zygosity"],
                "significance": h["significance"],
                "stars": h["stars"],
                "condition": h["conditions"][:48],
            }
            for h in r["hits"][: args.top]
        ]
        print(_table(rows, ["variant", "gene", "gt", "zygosity", "significance", "stars", "condition"]))
        print(f"  [{r['evidence']}]  {r['note']}")
        print(f"  stored under data/individuals/{args.name}/clinvar_screen.json (never under data/results)")
        return 0
    if args.action == "check":
        people = {p["name"]: p for p in ind.list_individuals()}
        if args.name not in people:
            print(f"{args.name}: not a local individual (genomeos individual list)")
            return 1
        chroms = args.chrom or [
            c for c in people[args.name]["chromosomes"] if (Path("data/reference") / f"{c}.fa").exists()
        ]
        if not chroms:
            print(f"{args.name}: no chromosome with both this person's rows and a local reference sequence")
            return 1
        worst = ""
        for chrom in chroms:
            r = ind.check(args.name, chrom)
            print(
                f"{chrom}: {r['variants']:,} variants ({r['snv']:,} SNVs, {r['phased']:,} phased); "
                f"hap1 applied {r['hap1']['applied']:,}, hap2 applied {r['hap2']['applied']:,}; "
                f"reference mismatches {r['reference_mismatches']:,} ({r['seconds']} s)"
            )
            worst = worst or ("" if r["verdict"].startswith("matches") else r["verdict"])
        print(f"  verdict: {worst or 'matches GRCh38 on every chromosome checked'}")
        print(f"  [{r['evidence']}]  stored under data/individuals/{args.name}/ (never under data/results)")
        return 0
    if args.action == "genes":
        from genomeos.genome import Annotation, IndexedGenome, default_gencode

        gff = default_gencode({args.chrom})
        fa = Path("data/reference") / f"{args.chrom}.fa"
        if not gff or not fa.exists():
            print(f"{args.chrom}: needs local models and sequence (genomeos data fetch --chrom {args.chrom})")
            return 1
        ann = Annotation.from_gff3(gff, {args.chrom})
        genome = IndexedGenome(str(fa))
        try:
            r = ind.gene_by_gene(args.name, args.chrom, ann, genome, limit=args.top)
        except FileNotFoundError as ex:
            print(ex)
            return 1
        finally:
            genome.close()
        if args.json:
            print(json.dumps(r, indent=2))
            return 0
        tt = r["totals"]
        print(
            f"{r['individual']} on {args.chrom}: {tt['genes_with_variants']:,} of {tt['genes']:,} coding "
            f"genes carry variants ({tt['variants_in_genes']:,} inside genes, {tt['coding_snvs']:,} coding "
            "SNVs)"
        )
        if r["by_consequence"]:
            print(
                "  coding SNVs by consequence: "
                + ", ".join(f"{k} {v}" for k, v in r["by_consequence"].items())
            )
        rows = [
            {
                "gene": g["gene"],
                "variants": f"{g['variants']:,}",
                "coding": g["coding_snvs"],
                "protein-changing": g["protein_changing"],
                "changes": "; ".join(
                    f"{x['hgvs_p'] or x['consequence']} ({x['genotype']})"
                    for x in g["changes"]
                    if x["consequence"] not in ("synonymous", "coding")
                )[:70],
            }
            for g in r["genes"]
        ]
        print(_table(rows, ["gene", "variants", "coding", "protein-changing", "changes"]))
        print(f"  [{r['evidence']}]")
        return 0
    for p in ind.list_individuals():
        tag = "built-in test human" if p["builtin"] else (p.get("note") or p.get("source_file", ""))
        n = f"{p['variants']:,} variants, " if p.get("variants") else ""
        print(f"{p['name']:12} {n}{len(p['chromosomes'])} chromosomes  {tag}")
        print(f"             [{p['evidence']}]")
        ck = ind.checks(p["name"])
        if ck:
            mism = sum(r["reference_mismatches"] for r in ck.values())
            bad = [c for c, r in ck.items() if not r["verdict"].startswith("matches")]
            state = "matches GRCh38" if not bad else f"assembly doubt on {', '.join(bad)}"
            print(f"             reference check: {len(ck)} chromosomes, {mism:,} mismatches, {state}")
    return 0


def _individual_predict(args: argparse.Namespace) -> int:
    """Feature d: what this person's regulatory variants do to one gene, predicted per variant."""
    from genomeos.genome import Annotation, IndexedGenome, default_gencode
    from genomeos.genome.individuals import vcf_path
    from genomeos.genome.regulation import regulation_of
    from genomeos.genome.regulatory import load_ccres
    from genomeos.predict import AlphaGenomeAdapter, status
    from genomeos.predict.individual_effects import predict_gene

    st = status()
    if not st["enabled"]:
        print(f"AlphaGenome predictions are disabled: {st['reason']}. To enable: {st['how']}")
        return 2
    vcf = vcf_path(args.name, args.chrom)
    if vcf is None:
        print(f"{args.name} has no rows on {args.chrom} (genomeos individual list)")
        return 1
    gff = default_gencode({args.chrom})
    fa = Path("data/reference") / f"{args.chrom}.fa"
    ccres = load_ccres(args.chrom)
    if not gff or not fa.exists() or not ccres:
        print(
            f"{args.chrom}: needs local models, sequence and elements "
            f"(genomeos data fetch --chrom {args.chrom})"
        )
        return 1
    ann = Annotation.from_gff3(gff, {args.chrom})
    genome = IndexedGenome(str(fa))
    try:
        reg = regulation_of(args.gene.upper(), args.chrom, ccres, ann, genome.lengths[args.chrom])
    except KeyError:
        print(f"{args.gene} not on {args.chrom}")
        return 1
    finally:
        genome.close()
    scorer = AlphaGenomeAdapter()._live_scorer(threshold=0.02)  # noqa: SLF001
    r = predict_gene(
        args.name,
        reg["gene"],
        args.chrom,
        reg,
        scorer,
        vcf,
        max_variants=args.max,
        progress=lambda m: print("  " + m, flush=True),
    )
    if args.json:
        print(json.dumps({"model": st["model"], **r}, indent=2))
        return 0
    print(
        f"{r['individual']} × {r['gene']}: {r['variants_in_elements']} variants in "
        f"{r['elements_considered']} elements that reach the gene; {r['variants_scored']} scored, "
        f"{r['variants_moving_gene']} move it by ≥ {r['min_effect_log2fc']} log2  [predicted, {st['model']}]"
    )
    h = r["haplotypes"]
    print(
        f"  haplotype sums: hap1 {h['hap1_log2_fold_change']:+.2f}, hap2 {h['hap2_log2_fold_change']:+.2f}"
        + (
            f", unphased heterozygous {', '.join(f'{x:+.2f}' for x in h['unphased_heterozygous_effects'])}"
            if h["unphased_heterozygous_effects"]
            else ""
        )
    )
    rows = [
        {
            "variant": f"{args.chrom}:{v['pos']:,} {v['ref']}>{v['alt']}",
            "gt": v["genotype"],
            "element": f"{v['element_class']} {v['element']}",
            "dist": f"{v['distance_to_tss']:,}",
            "log2FC": f"{v['effect']['log2_fold_change']:+.2f}",
            "tissue": (v["effect"]["tissue"] or "")[:32],
            "tracks": v["effect"]["tracks_moved"],
        }
        for v in r["variants"][: args.top]
    ]
    print(_table(rows, ["variant", "gt", "element", "dist", "log2FC", "tissue", "tracks"]))
    print(f"  [{r['evidence']}]  {h['note']}")
    return 0


def cmd_mouse(args: argparse.Namespace) -> int:
    """The second mammal: a mouse chromosome through the same code, its nodes held against ours."""
    from genomeos.genome.mouse import analyse
    from genomeos.results import save_result

    r = analyse(args.chrom, progress=lambda m: print("  " + m, flush=True))
    save_result(f"mouse_mm10_{args.chrom}", {k: v for k, v in r.items() if k != "mouse_domains"})
    c = r["node_comparison"]
    print(
        f"mouse {args.chrom} (mm10): {r['length']:,} bp, {r['genes']:,} genes "
        f"({r['coding_genes']:,} coding), {r['elements']:,} ENCODE elements, {r['domains']} nodes"
    )
    print(
        f"  {r['domains_with_coding_genes']} nodes with coding genes, median {r['domain_length_median']:,} bp"
    )
    frac = f"{c['fraction_conserved']:.0%}" if c["fraction_conserved"] is not None else "n/a"
    same = (
        f"{c['fraction_same_neighbourhood']:.0%}" if c["fraction_same_neighbourhood"] is not None else "n/a"
    )
    print(f"  against {r['human_nodes_indexed']:,} human nodes ({r['human_symbols_indexed']:,} symbols):")
    print(f"    {c['tested']} mouse nodes have ≥ 2 symbol matches; {c['unmapped']} have fewer")
    print(f"    {c['conserved']} land in one human node ({frac})")
    print(f"    {c['split_adjacent']} span adjacent human nodes, {c['split_scattered']} scattered")
    print(f"    same human neighbourhood (one node or adjacent nodes): {same}")
    rows = [
        {
            "mouse node": x["id"],
            "coding": x["coding_genes"],
            "mapped": x["mapped"],
            "verdict": x["verdict"],
            "human nodes": ", ".join(x["human_nodes"][:4]),
            "genes": ", ".join(x["genes"][:5]),
        }
        for x in sorted(r["node_rows"], key=lambda x: -x["mapped"])[: args.top]
    ]
    print(_table(rows, ["mouse node", "coding", "mapped", "verdict", "human nodes", "genes"]))
    for k, v in r["evidence"].items():
        print(f"  [{k}: {v}]")
    return 0


def cmd_lookup(args: argparse.Namespace) -> int:
    from genomeos.genome import Annotation, IndexedGenome, default_gencode
    from genomeos.genome.lookup import lookup, parse_variant

    chrom, pos, ref, alt = parse_variant(" ".join(args.variant))
    ann = genome = None
    gff = default_gencode({chrom})
    fa = Path(f"data/reference/{chrom}.fa.gz")
    if gff and fa.exists():
        ann = Annotation.from_gff3(gff, {chrom})
        genome = IndexedGenome(fa)
    r = lookup(chrom, pos, ref, alt, ann, genome, pathways=0 if args.no_pathways else 10)
    if genome:
        genome.close()
    v = r["vep"]
    damage = f"  SIFT {v['sift']}" if v["sift"] else ""
    damage += f"  PolyPhen {v['polyphen']}" if v["polyphen"] else ""
    print(r["variant"])
    print(f"  consequence: {v['consequence']} in {v['gene'] or '-'} {v['hgvsc']} {v['hgvsp']}{damage}")
    print(f"    [{v['evidence']}]")
    ids = f"  ({', '.join(v['dbsnp'])})" if v.get("dbsnp") else ""
    print("  known: " + "; ".join(r["known"]) + ids)
    if v.get("pubmed"):
        more = " …" if v["pubmed_count"] > 8 else ""
        print("  PubMed: " + ", ".join(str(p) for p in v["pubmed"][:8]) + more)
    for ind in r.get("individuals", []):
        gt = f" (genotype {ind['genotype']})" if ind["carries"] else ""
        has = "carries it" if ind["carries"] else "does not carry it"
        print(f"  {ind['individual']}: {has}{gt}  [{ind['evidence']}]")
    lt = r.get("local_trace")
    if lt:
        verdict = "agrees with VEP" if lt["agrees_with_vep"] else "differs from VEP (isoform or region)"
        names = f"{lt['hgvs_c'] or ''} {lt['hgvs_p'] or ''}"
        print(f"  local trace ({lt['transcript']}): {lt['consequence']} {names}")
        print(f"    {verdict}  [{lt['evidence']}]")
    p = r.get("protein")
    if p and "error" not in p:
        feats = "; ".join(
            f"{f['type']} {f['description']} ({f['start']}-{f['end']})" for f in p["features_at_residue"]
        )
        iso = "" if p["reference_residue_matches"] else " (reference residue differs: isoform?)"
        print(f"  protein {p['accession']} {p['name']} ({p['length']} aa), residue {p['residue']}{iso}:")
        af = ", AlphaFold" if p["alphafold"] else ""
        nx = p["structures_experimental"]
        print(f"    {feats or 'no annotated feature'}; structures: {nx} experimental{af}")
        print(f"    [{p['evidence']}]")
    elif p:
        print(f"  protein: {p['error']}")
    pw = r.get("pathways")
    if pw:
        w = pw["most_affected"]
        worst = f"; most affected {w['name']} ({w['fraction_lost']:.0%})" if w else ""
        print(
            f"  pathways: {pw['reactions_lost']} reactions lost over {pw['checked']} pathways with the "
            f"protein {pw['treated_as']}{worst}"
        )
        print(f"    [{pw['evidence']}]")
    if args.json:
        print(json.dumps(r, indent=1))
    return 0


def cmd_predict(args: argparse.Namespace) -> int:
    """Optional AlphaGenome feature: loaded always, enabled only with the package and the key."""
    import re

    from genomeos.predict import AlphaGenomeAdapter, status

    st = status()
    if args.status or (not args.variant and not args.element):
        why = f" ({st['reason']})" if st["reason"] else ""
        print(f"{st['name']}: {'enabled' if st['enabled'] else 'disabled'}{why}")
        for f in st["features"]:
            print(f"  {f['id']}. {f['name']:18} {f['state']:8} {f['what']}")
        if not st["enabled"]:
            print(f"  to enable: {st['how']}")
        print(f"  {st['licence']}")
        return 0
    if not st["enabled"]:
        print(f"AlphaGenome predictions are disabled: {st['reason']}. To enable: {st['how']}")
        return 2
    if args.element:
        return _predict_element(args, st)
    m = re.match(r"^(chr\w+):(\d+)\s+([ACGTacgt]+)>([ACGTacgt]+)$", " ".join(args.variant).strip())
    if not m:
        print("expected a variant like chr21:25897620 C>T (1-based, plus strand)")
        return 2
    chrom, pos, ref, alt = m.group(1), int(m.group(2)), m.group(3).upper(), m.group(4).upper()
    adapter = AlphaGenomeAdapter()
    effects = adapter.predict(chrom, pos, ref, alt, threshold=args.min)
    effects.sort(key=lambda e: -abs(e.log2_fold_change))
    scan = adapter.last_scan
    if args.json:
        print(
            json.dumps(
                {
                    "variant": f"{chrom}:{pos} {ref}>{alt}",
                    "model": st["model"],
                    "evidence": "predicted",
                    "threshold_log2fc": args.min,
                    "scanned": scan,
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
                },
                indent=2,
            )
        )
        return 0
    print(f"{chrom}:{pos} {ref}>{alt}: {len(effects)} predicted tissue effects [predicted, {st['model']}]")
    if scan:
        print(
            f"  scanned {scan['genes']} genes in the 1 Mb window across {scan['tracks']} RNA-seq tracks; "
            f"largest |log2 fold change| {scan['max_abs_log2fc']:.3f} (threshold {args.min})"
        )
    for e in effects[: args.top]:
        d = getattr(e.direction, "value", str(e.direction))
        print(
            f"  {e.gene:12} {e.tissue[:44]:44} log2FC {e.log2_fold_change:+.2f}  {d:9} "
            f"confidence {e.confidence:.2f}"
        )
    if len(effects) > args.top:
        print(f"  … {len(effects) - args.top} more (--top N or --json)")
    return 0


def _predict_element(args: argparse.Namespace, st: dict) -> int:
    """Feature b: delete an ENCODE element (or any region) and read which gene moves."""
    import re

    from genomeos.predict import AlphaGenomeAdapter
    from genomeos.predict.enhancer_target import Context

    m = re.match(r"^(chr\w+):([\d,]+)-([\d,]+)$", args.element.strip())
    if not m:
        print("expected a region like chr21:42721378-42721727 (0-based start, half-open) or an element id")
        return 2
    chrom, start, end = m.group(1), int(m.group(2).replace(",", "")), int(m.group(3).replace(",", ""))
    try:
        ctx = Context(chrom)
    except FileNotFoundError as ex:
        print(ex)
        return 1
    try:
        r = ctx.score_region(
            AlphaGenomeAdapter()._live_scorer(threshold=0.0), start, end, min_effect=args.min
        )  # noqa: SLF001
    finally:
        ctx.close()
    if args.json:
        print(json.dumps({"model": st["model"], "evidence": "predicted", **r}, indent=2))
        return 0
    print(f"{r['id']} {chrom}:{r['start']:,}-{r['end']:,} ({r['length']} bp) deleted in its 1 Mb window")
    print(
        f"  {r['genes_in_window']} genes read across {r['tracks']} RNA-seq tracks  [predicted, {st['model']}]"
    )
    if r.get("inferred"):
        print(
            f"  inferred target (nearest TSS in domain): {r['inferred']['gene']} "
            f"at {r['inferred']['distance']:,} bp"
        )
    p = r["predicted"]
    if p:
        print(
            f"  predicted target: {p['gene']} ({p['action']}, log2FC {p['log2_fold_change']:+.2f} "
            f"in {p['tissue']}, {p['strength']}, confidence {p['confidence']})"
        )
    else:
        print(f"  predicted target: none (no gene moves by {args.min} log2)")
    pc = r.get("predicted_coding")
    if pc and (not p or pc["gene"] != p["gene"]):
        print(
            f"  strongest coding gene: {pc['gene']} ({pc['action']}, {pc['log2_fold_change']:+.2f} "
            f"in {pc['tissue']})"
        )
    if r.get("verdict_coding"):
        print(f"  verdict: {r['verdict_coding']}")
    print("  genes that move most:")
    for g in r["top_genes"]:
        print(f"    {g['gene']:16} drop {g['max_drop_log2fc']:+.3f} in {g['max_drop_tissue'][:40]}")
    return 0


def cmd_ptm(args: argparse.Namespace) -> int:
    """Post-translational state: modifiable sites and their writers, per protein or genome-wide."""
    from genomeos.molecules.ptm import CACHE, INDEX, build_index, load_index, sites, summary
    from genomeos.results import save_result

    if args.protein:
        p = CACHE / f"{args.protein.upper()}.json"
        if not p.exists():
            print(f"{args.protein}: no compiled definition (genomeos protein {args.protein} --compile)")
            return 1
        ss = sites(json.loads(p.read_text()))
        by = {}
        for s in ss:
            by[s["class"]] = by.get(s["class"], 0) + 1
        print(
            f"{args.protein.upper()}: {len(ss)} modifiable sites  "
            + ", ".join(f"{k} {v}" for k, v in by.items())
        )
        rows = [
            {
                "site": f"{s['start']}" if s["start"] == s["end"] else f"{s['start']}-{s['end']}",
                "class": s["class"],
                "modification": s["description"][:44],
                "written by": ", ".join(s["writers"]) or "",
            }
            for s in ss[: args.top]
        ]
        print(_table(rows, ["site", "class", "modification", "written by"]))
        print("  [curated: UniProt features; a site that can be modified, not a measurement that it is]")
        return 0
    idx = load_index(INDEX) if not args.rebuild else None
    if idx is None:
        idx = build_index()
    if args.writer:
        subs = idx["writers"].get(args.writer.upper(), {})
        print(
            f"{args.writer.upper()}: {len(subs)} substrates with a UniProt-recorded site  [curated: UniProt]"
        )
        for g, n in sorted(subs.items(), key=lambda kv: (-kv[1], kv[0]))[: args.top]:
            print(f"  {g:12} {n} site{'s' if n != 1 else ''}")
        return 0
    s = summary(idx, top=args.top)
    save_result("ptm_genome_wide", s)
    print(
        f"{s['proteins_with_sites']:,} of {s['proteins']:,} compiled proteins carry "
        f"{s['sites']:,} modifiable sites; {s['writers']} named writers, "
        f"{s['writer_edges']:,} writer→substrate edges  [{s['evidence']}]"
    )
    print("  by class: " + ", ".join(f"{k} {v:,}" for k, v in s["by_class"].items()))
    print(_table(s["top_writers"], ["writer", "substrates", "sites"]))
    return 0


def cmd_graph(args: argparse.Namespace) -> int:
    from genomeos.molecules.graph import build, summarise
    from genomeos.results import save_result

    g = build(min_score=args.min_score)
    if args.gene:
        n = g.neighbourhood(args.gene.upper(), args.max_nodes)
        if not n["nodes"]:
            print(f"{args.gene} is not in the local knowledge graph; compile it first:")
            print(f"  genomeos protein {args.gene} --compile")
            return 1
        c = n["centre"]
        print(f"{c}: {len(n['nodes']) - 1} neighbours, {len(n['edges'])} edges in the neighbourhood")
        rows = []
        for e in sorted(n["edges"], key=lambda e: -e["confidence"]):
            if c not in (e["a"], e["b"]):
                continue
            other = e["b"] if e["a"] == c else e["a"]
            node = g.nodes[other]
            detail = f"score {e['score']}" if "score" in e else ""
            if "ntpm" in e:
                detail = f"{e['ntpm']:.0f} nTPM"
            rows.append(
                {
                    "relation": e["rel"],
                    "node": node.get("name") or other,
                    "kind": node["kind"],
                    "evidence": e["evidence"],
                    "conf": e["confidence"],
                    "detail": detail,
                }
            )
        print(_table(rows[: args.top], ["relation", "node", "kind", "evidence", "conf", "detail"]))
        return 0
    s = summarise(g)
    print(f"knowledge graph: {s['nodes']:,} nodes {s['node_kinds']}, {s['edges']:,} edges {s['edge_kinds']}")
    print(
        f"  compiled proteins {s['compiled_proteins']}, physical associations "
        f"{s['physical_associations']:,}, isolated {s['isolated_proteins']}, components {s['components']} "
        f"(largest {s['largest_component']})"
    )
    print("  hubs: " + ", ".join(f"{h['gene']} {h['associations']}" for h in s["hubs"][:8]))
    print("  biggest pathways: " + "; ".join(f"{n} ({c})" for n, c in s["biggest_pathways"][:6]))
    print(f"  [{s['evidence']}]")
    if args.save:
        save_result(f"graph_{args.save}", s)
        print(f"  saved data/results/graph_{args.save}.json")
    return 0


def _rna_isoforms(args: argparse.Namespace) -> int:
    """Isoform-level expression: which transcript each tissue actually makes (GTEx transcript TPM)."""
    from genomeos.molecules.rna import gtex_isoforms

    sym = args.symbol.upper()
    transcripts = []
    p = Path("data/knowledge/proteins") / f"{sym}.json"
    if p.exists():
        d = json.loads(p.read_text())
        transcripts = ((d["sections"].get("genomic_origin") or {}).get("items") or {}).get(
            "transcripts"
        ) or []
    r = gtex_isoforms(sym, transcripts=transcripts)
    if r.get("error"):
        print(f"{sym}: {r['error']}")
        return 1
    print(
        f"{sym}: {r['transcripts_measured']} transcripts measured in {r['tissues_measured']} tissues; "
        f"canonical {r['canonical'] or '?'} is the dominant isoform in {r['canonical_dominant_in']} tissues, "
        f"another in {len(r['tissues_where_another_isoform_dominates'])}  [{r['evidence']}]"
    )
    rows = [
        {
            "transcript": tid,
            "name": v["name"] or "",
            "canonical": "yes" if v["canonical"] else "",
            "aa": v["protein_length"] or "",
            "dominant in": v["tissues_dominant"],
            "median TPM": f"{v['median_tpm']:.2f}",
            "max TPM": f"{v['max_tpm']:.2f}",
        }
        for tid, v in list(r["isoforms"].items())[: args.top]
    ]
    print(_table(rows, ["transcript", "name", "canonical", "aa", "dominant in", "median TPM", "max TPM"]))
    other = r["tissues_where_another_isoform_dominates"]
    if other:
        print(
            "  tissues where a non-canonical isoform dominates: "
            + "; ".join(
                f"{t} → {d['name'] or d['transcript']} ({d['share']:.0%})" for t, d in list(other.items())[:8]
            )
        )
    return 0


def cmd_rna(args: argparse.Namespace) -> int:
    if getattr(args, "isoforms", False):
        return _rna_isoforms(args)
    from genomeos.genome import Annotation, IndexedGenome, default_gencode
    from genomeos.molecules.rna import rna_report

    ann = genome = None
    gff = default_gencode({args.chrom}) if args.chrom else None
    if gff:
        ann = Annotation.from_gff3(gff, {args.chrom})
        if Path(args.genome).exists():
            genome = IndexedGenome(args.genome)
    r = rna_report(args.symbol, ann, genome, expression=not args.no_expression)
    if genome:
        genome.close()
    tx = r.get("transcripts")
    if tx:
        print(
            f"{tx['gene']} ({tx['gene_type']}): {tx['count']} transcripts, {tx['coding_isoforms']} coding; "
            f"biotypes {tx['by_biotype']}  [{tx['evidence']}]"
        )
        print(
            _table(
                [
                    {
                        "transcript": x["name"],
                        "biotype": x["biotype"],
                        "exons": x["exons"],
                        "spliced nt": x.get("spliced_nt", "-"),
                        "CDS nt": x.get("cds_nt", "-"),
                        "aa": x.get("protein_aa", "-"),
                        "tags": ", ".join(x["tags"]),
                    }
                    for x in tx["transcripts"][: args.top]
                ],
                ["transcript", "biotype", "exons", "spliced nt", "CDS nt", "aa", "tags"],
            )
        )
    elif ann is not None:
        print(f"{args.symbol} not on {args.chrom}")
    e = r.get("expression")
    if e and e.get("tissues"):
        print(
            f"expression (GTEx): {e['pattern']}; median {e['median_tpm']} TPM over "
            f"{e['tissues_measured']} tissues, {e['tissues_expressed']} with ≥1 TPM  "
            f"[{e['evidence']}] conf={e['confidence']}"
        )
        print("  top: " + ", ".join(f"{t_} {v:.0f}" for t_, v in e["top"]))
    elif e:
        print(f"expression: {e.get('error', 'unavailable')}")
    return 0


def cmd_repeats(args: argparse.Namespace) -> int:
    from genomeos.genome.repeats import fetch_repeats, load_repeats, save_repeats, summarise
    from genomeos.results import save_result

    reps = load_repeats(args.chrom) if not args.refresh else []
    if not reps:
        print(f"fetching RepeatMasker rows for {args.chrom} from UCSC…", flush=True)
        reps = fetch_repeats(args.chrom)
        save_repeats(args.chrom, reps)
    length = None
    if Path(args.genome).exists():
        from genomeos.genome import IndexedGenome

        g = IndexedGenome(args.genome)
        length = g.lengths.get(args.chrom)
        g.close()
    s = summarise(args.chrom, reps, length)
    save_result(f"rmsk_{args.chrom}", s)
    print(
        f"{args.chrom}: {s['copies']:,} repeat copies, {s['repeat_bp'] / 1e6:.1f} Mb"
        + (f" ({s['repeat_fraction']:.1%} of the chromosome)" if s["repeat_fraction"] else "")
        + f"  [curated: {s['evidence']}]"
    )
    print(
        _table(
            [
                {
                    "class": k,
                    "copies": f"{v['copies']:,}",
                    "Mb": f"{v['bp'] / 1e6:.2f}",
                    "mean divergence": f"{v['mean_divergence']:.0%}",
                }
                for k, v in s["by_class"].items()
            ],
            ["class", "copies", "Mb", "mean divergence"],
        )
    )
    print(
        "  top families: "
        + ", ".join(f"{k} {v / 1e6:.1f} Mb" for k, v in list(s["top_families"].items())[:8])
    )
    print(f"  saved data/results/rmsk_{args.chrom}.bed.gz and rmsk_{args.chrom}.json")
    return 0


def _predicted_cell(e: dict) -> str:
    p = e.get("predicted")
    if p:
        return f"{p['gene']} {p['log2_fold_change']:+.2f} {p['tissue'][:18]}"
    return "no effect" if "predicted" in e else ""


def cmd_regulation(args: argparse.Namespace) -> int:
    from genomeos.genome import Annotation, Genome, default_gencode
    from genomeos.genome.regulation import regulation_of
    from genomeos.genome.regulatory import load_ccres

    ccres = load_ccres(args.chrom)
    if not ccres:
        print(f"no ENCODE elements for {args.chrom}")
        return 1
    genome = Genome.from_fasta(args.genome)
    ann = Annotation.from_gff3(args.gff3 or default_gencode({args.chrom}), {args.chrom})
    try:
        r = regulation_of(args.gene, args.chrom, ccres, ann, genome.chromosomes[args.chrom].length)
    except KeyError:
        print(f"{args.gene} not on {args.chrom}")
        return 1
    d = r["domain"]
    print(f"{r['gene']} TSS {args.chrom}:{r['tss']:,} ({r['strand']})")
    if d:
        print(
            f"  node {d['id']}: {d['start']:,}-{d['end']:,} ({d['length']:,} bp), {d['coding_genes']} coding "
            f"genes, {d['enhancers']} enhancer-like elements  [{d['evidence'][:40]}…] conf={d['confidence']}"
        )
        print(f"  genes sharing the node: {', '.join(d['genes'][:10])}{'…' if len(d['genes']) > 10 else ''}")
    print(
        f"  promoter elements: {len(r['promoters'])}  "
        + "; ".join(
            f"{p['id']} ({p['class']}, {p['distance']} bp, conf {p['confidence']})" for p in r["promoters"]
        )
    )
    print(
        f"  enhancers that can reach it: {r['enhancers_in_domain']} in the node, "
        f"{r['enhancers_nearest_to_this_gene']} nearest to this gene, "
        f"{r['enhancers_inside_gene']} inside the gene"
    )
    if r.get("enhancers_predicted"):
        print(
            f"  predicted (AlphaGenome, deletion): {r['enhancers_predicted']} scored, "
            f"{r['enhancers_predicted_this_gene']} name this gene"
        )
    print(f"  insulators bounding the node: {len(r['insulators_bounding'])}")
    print(
        _table(
            [
                {
                    "element": e["id"],
                    "class": e["class"],
                    "position": f"{e['start']:,}",
                    "distance": f"{e['distance']:,}",
                    "where": "intragenic" if e["intragenic"] else "outside",
                    "basis": e["basis"],
                    "conf": e["confidence"],
                    "predicted": _predicted_cell(e),
                }
                for e in r["enhancers"][: args.top]
            ],
            ["element", "class", "position", "distance", "where", "basis", "conf", "predicted"],
        )
    )
    for k, v in r["evidence"].items():
        print(f"  [{k}: {v}]")
    return 0


def cmd_domains(args: argparse.Namespace) -> int:
    from genomeos.genome import Annotation, Genome, default_gencode
    from genomeos.genome.domains import infer_domains, summarise
    from genomeos.genome.regulatory import load_ccres
    from genomeos.results import save_result

    ccres = load_ccres(args.chrom)
    if not ccres:
        print(
            f"no ENCODE elements for {args.chrom}; run `genomeos data distil --only encode_ccres_chr21` "
            "or the genome-wide job"
        )
        return 1
    genome = Genome.from_fasta(args.genome)
    length = genome.chromosomes[args.chrom].length
    ann = Annotation.from_gff3(args.gff3 or default_gencode({args.chrom}), {args.chrom})
    doms = infer_domains(args.chrom, length, ccres, ann)
    s = summarise(doms)
    if getattr(args, "hic", None):
        from genomeos.genome import hic
        from genomeos.predict.contact_maps import compare, random_control

        st = hic.status()
        if not st["enabled"]:
            print(f"4DN boundary downloads are disabled: {st['how']}")
            return 2
        try:
            man = hic.fetch_boundaries(args.hic, progress=lambda m: print("  " + m, flush=True))
        except (PermissionError, LookupError, OSError) as ex:
            print(ex)
            return 1
        measured = hic.load_boundaries(args.hic, args.chrom)
        inferred = sorted({d.start for d in doms} | {d.end for d in doms})
        inferred = [b for b in inferred if 0 < b < length]
        cmp = compare(inferred, measured)
        cmp["random_control"] = random_control(inferred, measured, length)
        save_result(
            f"domains_{args.chrom}_hic_{args.hic.replace(' ', '_').replace('/', '_')}",
            {
                "chrom": args.chrom,
                "biosource": args.hic,
                "file": man.get("accession"),
                "comparison": cmp,
                "evidence": {
                    "inferred": "inferred: CTCF-only ENCODE elements as boundaries",
                    "measured": st["evidence"],
                },
            },
        )
        print(
            f"{args.chrom}: {cmp['inferred']} CTCF-only boundaries against {cmp['predicted']} measured "
            f"{args.hic} boundaries ({man.get('accession')}, {man.get('experiment_type')})"
        )
        print(
            f"  {cmp['inferred_on_a_predicted_boundary']} inferred boundaries sit within "
            f"{cmp['tolerance'] // 1000} kb of a measured one ({cmp['fraction_inferred_supported']:.0%}, "
            f"against {cmp['random_control']:.0%} at random); {cmp['predicted_on_an_inferred_boundary']} "
            f"measured boundaries have an inferred one ({cmp['fraction_predicted_matched']:.0%}); "
            f"median distance {cmp['median_distance_inferred_to_predicted']:,} bp"
        )
        return 0
    if getattr(args, "contact_map", False):
        from genomeos.predict import status
        from genomeos.predict.contact_maps import Insulation, compare
        from genomeos.predict.splice_sites import client_factory

        st = status()
        if not st["enabled"]:
            print(f"AlphaGenome contact maps are disabled: {st['reason']}. To enable: {st['how']}")
            return 2
        ins = Insulation(args.chrom, length).load(
            client_factory, progress=lambda m: print("  " + m, flush=True)
        )
        predicted = ins.boundaries()
        inferred = sorted({d.start for d in doms} | {d.end for d in doms})
        inferred = [b for b in inferred if 0 < b < length]
        from genomeos.predict.contact_maps import random_control

        cmp = compare(inferred, [b["pos"] for b in predicted])
        cmp["random_control"] = random_control(inferred, [b["pos"] for b in predicted], length)
        save_result(
            f"domains_{args.chrom}_contact",
            {
                "chrom": args.chrom,
                "windows": len(ins.windows),
                "tracks_per_window": ins.windows[0]["tracks"] if ins.windows else 0,
                "predicted_boundaries": predicted,
                "comparison": cmp,
                "evidence": {
                    "inferred": "inferred: CTCF-only ENCODE elements as boundaries",
                    "predicted": (
                        f"predicted: {st['model']} contact maps, insulation minima at 2 kb, 28 cell types"
                    ),
                },
            },
        )
        print(
            f"{args.chrom}: {cmp['inferred']} CTCF-only boundaries against {cmp['predicted']} insulation "
            f"minima from the predicted contact map ({len(ins.windows)} windows, "
            f"{ins.windows[0]['tracks']} cell types):"
        )
        print(
            f"  {cmp['inferred_on_a_predicted_boundary']} inferred boundaries sit within "
            f"{cmp['tolerance'] // 1000} kb of a predicted one ({cmp['fraction_inferred_supported']:.0%}, "
            f"against {cmp['random_control']:.0%} for boundaries placed at random); "
            f"{cmp['predicted_on_an_inferred_boundary']} predicted minima have an inferred boundary "
            f"({cmp['fraction_predicted_matched']:.0%}); median distance "
            f"{cmp['median_distance_inferred_to_predicted']:,} bp"
        )
        print(
            "  [predicted evidence, capped at 0.7; where the two agree the node's edge has two kinds "
            "of evidence]"
        )
        return 0
    print(
        f"{args.chrom}: {s['domains']} domains from CTCF-only boundaries; "
        f"median {s['size_median']:,} bp, max {s['size_max']:,} bp"
    )
    print(
        f"  coding genes per domain: median {s['coding_genes_per_domain_median']}, "
        f"max {s['largest_gene_count']}; {s['domains_without_coding_genes']} domains without coding genes; "
        f"enhancers per domain median {s['enhancers_per_domain_median']}"
    )
    print(f"  [{s['evidence']}] conf={s['confidence']}")
    print(
        _table(
            [
                {
                    "domain": d.id,
                    "start": f"{d.start:,}",
                    "size": f"{d.length:,}",
                    "coding": d.coding_genes,
                    "promoters": d.promoters,
                    "enhancers": d.enhancers,
                    "genes": (", ".join(d.genes[:6]) + ("…" if len(d.genes) > 6 else "")),
                }
                for d in sorted(doms, key=lambda d: -d.coding_genes)[: args.top]
            ],
            ["domain", "start", "size", "coding", "promoters", "enhancers", "genes"],
        )
    )
    save_result(f"domains_{args.chrom}", {**s, "domains": [d.to_dict(compact=True) for d in doms]})
    print(f"  saved data/results/domains_{args.chrom}.json")
    return 0


def cmd_libs(args: argparse.Namespace) -> int:
    if getattr(args, "proteome", False):
        from genomeos.lib.proteome import summary

        s = summary()
        print(f"packaged proteome: {s['proteins']:,} human proteins (genomeos/lib/data/proteome.json.gz)")
        for k in (
            "with_function",
            "with_pathways",
            "with_partners",
            "with_experimental_structure",
            "with_disease",
        ):
            print(f"  {k.replace('_', ' '):<28}{s[k]:>7,}  {s[k] / s['proteins']:.1%}")
        print(f"  {'symbol mismatch (marked)':<28}{s['symbol_mismatch']:>7,}")
        for k, v in s["evidence"].items():
            print(f"  [{k}: {v}]")
        return 0
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


def cmd_therapeutic(args: argparse.Namespace) -> int:
    """Cancer targeting and therapeutic mechanism reasoning for one tumour."""
    from genomeos.therapeutics import DISCLAIMER, analyse_vcf, text_report, write_outputs
    from genomeos.therapeutics.report import target_specification_text

    if args.install_hla_models:
        from genomeos.therapeutics.binding import fetch_mhcflurry_models

        print(fetch_mhcflurry_models(log=sys.stdout))
        return 0

    hla = [a for spec in (args.hla or []) for a in spec.split(",") if a.strip()]
    a = analyse_vcf(
        args.tumour,
        normal_vcf=args.normal,
        hla=hla,
        rna=args.rna,
        chroms=set(args.chrom) if args.chrom else None,
        top_genes=args.top,
        purity=args.purity,
        cnv=args.cnv,
        cohort=args.cohort or "",
        net=not args.offline,
        indirect=not args.no_indirect,
        hla_predictor=args.hla_predictor,
        log=sys.stdout,
    )
    if args.report:
        print(text_report(a, detail=args.detail))
    else:
        print(
            f"{a['sample']}: {a['variants_total']} variants, {a['somatic_candidates']} somatic "
            f"candidates, {len(a['candidates'])} target candidates "
            f"(data level {a['data_level']['level_reached']}/9)"
        )
        rows = []
        for c in a["candidates"]:
            best = c.best_mechanism
            rows.append(
                {
                    "gene": c.gene,
                    "class": c.target_class.replace("_", " "),
                    "score": "-" if c.scores.overall is None else f"{c.scores.overall:.2f}",
                    "access": _fmt(c.scores.value("surface_accessibility")),
                    "select": _fmt(c.scores.value("tumour_selectivity")),
                    "safety": _fmt(c.scores.value("normal_tissue_safety")),
                    "intern": _fmt(c.scores.value("internalisation")),
                    "best mechanism": f"{best.mechanism} {best.compatibility:.0%}" if best else "none",
                }
            )
        print(
            _table(
                rows,
                ["gene", "class", "score", "access", "select", "safety", "intern", "best mechanism"],
            )
        )
        print(f"\n{DISCLAIMER}")
    if args.spec:
        for c in a["candidates"]:
            if c.gene.upper() == args.spec.upper():
                print(target_specification_text(c, a.get("provider_bundle")))
                break
        else:
            print(f"{args.spec} is not among the candidates")
    if args.out:
        files = write_outputs(a, args.out)
        print("\nwritten:")
        for path in files.values():
            print(f"  {path}")
    return 0


def _fmt(v: float | None) -> str:
    return "-" if v is None else f"{v:.2f}"


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
    p.add_argument("--bai", help="index of an indexed BAM given as `reads` (URL or path): read it by ranges")
    p.add_argument("--chrom", nargs="*", help="chromosome ends to read (default: autosomes and chrX)")
    p.add_argument("--window", type=int, default=10_000, help="bases read at each chromosome end")
    p.add_argument("--sample-mb", type=float, default=100.0, help="megabytes of the unmapped tail to sample")
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

    p = sub.add_parser("grow", help="grow an organism from one cell (BioLang v0.3 organism program)")
    p.add_argument("module", help="a .bio file with an organism block")
    p.add_argument("--until", default="800", help="organism time: 800, '14 h', '20 yr' (default 800 min)")
    p.add_argument("--seed", type=int, default=None, help="seed for timer spread (default: the program's)")
    p.add_argument("--means", action="store_true", help="run every timer at its mean, no spread")
    p.add_argument("--depth", type=int, default=2, help="lineage tree depth to print (0 = none)")
    p.add_argument("--max-cells", type=int, default=200_000)
    p.add_argument("--compare", action="store_true", help="score against the organism's reference lineage")
    p.add_argument("--json", help="write summary, asserts, uncertainty and diff to this file")
    p.add_argument(
        "--experiments",
        nargs="?",
        const="",
        default=None,
        metavar="NAMES",
        help="run the program's experiment blocks (all, or a comma-separated list) against the wild type",
    )
    p.add_argument(
        "--design",
        nargs="?",
        const="",
        default=None,
        metavar="NAMES",
        help="run the program's design blocks with BioForge (all, or a comma-separated list)",
    )
    p.set_defaults(fn=cmd_grow)

    p = sub.add_parser("organism", help="minimal organism: C. elegans early lineage")
    p.add_argument("--minutes", type=float, default=150.0)
    p.add_argument("--depth", type=int, default=2)
    p.set_defaults(fn=cmd_organism)

    data = sub.add_parser("data", help="storage economics: status, distil, clean, results").add_subparsers(
        dest="data_cmd", required=True
    )
    data.add_parser("status", help="disk use and what can be distilled")
    p = data.add_parser(
        "fetch", help="bring one more chromosome to full footing: sequence, models, elements, repeats"
    )
    p.add_argument("--chrom", nargs="+", required=True)
    p.add_argument("--no-elements", action="store_true")
    p.add_argument("--no-repeats", action="store_true")
    p.add_argument("--individual", action="store_true", help="also the test human's (HG002) variants")
    p.add_argument("--analyse", action="store_true", help="then classify the UNKNOWN space and infer domains")
    p = data.add_parser("distil", help="turn raw downloads into result summaries")
    p.add_argument("--only", nargs="*")
    p.add_argument("--force", action="store_true", help="recompute existing summaries")
    p = data.add_parser("clean", help="delete raw inputs that already have a summary")
    p.add_argument("--yes", action="store_true", help="actually delete (default is a dry run)")
    data.add_parser("results", help="list result summaries")
    for sp in data.choices.values():
        sp.set_defaults(fn=cmd_data)

    p = sub.add_parser("calibrate", help="fit an ageing parameter to a measurement (rules from data)")
    p.add_argument("--cell-type", default="hematopoietic_stem", choices=sorted(CELL_TYPES))
    p.add_argument("--age", type=float, required=True)
    p.add_argument("--measured", type=float, required=True, help="mean telomere length, bp (Southern scale)")
    p.add_argument("--source", required=True, help="where the measurement comes from")
    p.set_defaults(fn=cmd_calibrate)

    p = sub.add_parser("methylation", help="clock CpG methylation from bedMethyl files or URLs (streamed)")
    p.add_argument("sources", nargs="+", help="bedMethyl paths or http(s) URLs; haplotype files are pooled")
    p.add_argument("--min-coverage", type=int, default=5)
    p.add_argument("--save", metavar="NAME")
    p.set_defaults(fn=cmd_methylation)

    sig = sub.add_parser(
        "signals", help="the genome's delimiters: learn PWMs from annotation, scan raw DNA"
    ).add_subparsers(dest="signals_cmd", required=True)
    p = sig.add_parser("learn", help="learn splice/start signals from a chromosome's annotation")
    p.add_argument("--chrom", default="chr21")
    p.add_argument("--genome", default="data/reference/chr21.fa.gz")
    p.add_argument("--gff3")
    p.add_argument("-o", "--output")
    p.set_defaults(fn=cmd_signals)
    p = sig.add_parser("scan", help="score a locus with learned signals")
    p.add_argument("locus")
    p.add_argument("--signals", default="data/results/signals_chr21.json")
    p.add_argument("--genome", default="data/reference/chr21.fa.gz")
    p.add_argument("--min-relative", type=float, default=0.8)
    p.add_argument("--limit", type=int, default=40)
    p.set_defaults(fn=cmd_signals)

    p = sub.add_parser(
        "anatomy",
        help="count the blocks and elements of a chromosome; --compare shows saved genomes side by side",
    )
    p.add_argument("fasta", nargs="?")
    p.add_argument("--chrom")
    p.add_argument("--gff3")
    p.add_argument("--name")
    p.add_argument("--save", metavar="NAME")
    p.add_argument("--compare", action="store_true")
    p.set_defaults(fn=cmd_anatomy)

    p = sub.add_parser("design", help="parts list and base budget for a cell type built from the libraries")
    p.add_argument("cell_type")
    p.add_argument("--save", action="store_true")
    p.set_defaults(fn=cmd_design)

    can = sub.add_parser(
        "cancer", help="healthy vs tumour comparison and cBioPortal-distilled cancer knowledge"
    ).add_subparsers(dest="cancer_cmd", required=True)
    p = can.add_parser("distil", help="distil driver-gene frequencies from a cBioPortal study")
    p.add_argument("--study", default="msk_impact_2017")
    p = can.add_parser(
        "expression",
        help="distil a study's tumour-versus-normal expression per gene (cohort reference)",
    )
    p.add_argument("--study", required=True, help="cBioPortal study id")
    p.add_argument("--genes", nargs="*", help="gene symbols; the driver panel by default")
    p.add_argument("--top", type=int, default=20, help="rows to print")
    p = can.add_parser("genes", help="most frequently mutated driver genes")
    p.add_argument("--top", type=int, default=25)
    p = can.add_parser("gene", help="one gene: frequency, hotspots, cancer types")
    p.add_argument("symbol")
    p = can.add_parser("types", help="cancer types in the distilled study")
    p.add_argument("--top", type=int, default=30)
    p = can.add_parser(
        "tumour", help="tumour DNA alone: VEP annotation, germline estimate, drivers, peptides, pathways"
    )
    p.add_argument("--vcf", required=True)
    p.add_argument("--chrom", nargs="*")
    p.add_argument(
        "--deep", type=int, default=5, help="top coding variants to trace to peptides and pathways"
    )
    p.add_argument("--no-pathways", action="store_true")
    p.add_argument(
        "--therapeutic",
        action="store_true",
        help="also run the therapeutic target pipeline (genomeos therapeutic has the full report)",
    )
    p.add_argument("--packet", help="write the AI-agent task packet (JSON)")
    p = can.add_parser("compare", help="somatic variants of a tumour vs the same person's normal sample")
    p.add_argument("--normal", required=True)
    p.add_argument("--tumour", required=True)
    p.add_argument("--genome", default="data/reference/chr21.fa.gz")
    p.add_argument("--gff3")
    p.add_argument("--chrom", nargs="*")
    p.add_argument("--packet", help="write the AI-agent task packet (JSON)")
    for sp in can.choices.values():
        sp.set_defaults(fn=cmd_cancer)

    p = sub.add_parser(
        "unknown", help="investigate UNKNOWN blocks: classify the space between genes, largest first"
    )
    p.add_argument("--genome", default="data/reference/chr21.fa.gz")
    p.add_argument("--chrom")
    p.add_argument("--gff3")
    p.add_argument("--patterns", help="TSV of named motifs (name, class, regex, evidence, confidence, note)")
    p.add_argument("--min-size", type=int, default=1000)
    p.add_argument("--max-blocks", type=int)
    p.add_argument("--save")
    p.set_defaults(fn=cmd_unknown)

    p = sub.add_parser(
        "budget", help="the 98%%: every UNKNOWN block with constraint attached and a best guess"
    )
    p.add_argument("--chrom", help="one chromosome (reads Zoonomia phyloP over its blocks); omit: the genome")
    p.add_argument(
        "--threshold", type=float, default=2.27, help="phyloP for a constrained base (Zoonomia 5%% FDR)"
    )
    p.add_argument("--no-phylop", action="store_true")
    p.add_argument("--no-elements", action="store_true")
    p.add_argument("--cached", action="store_true", help="show the saved budget instead of recomputing")
    p.add_argument("--top", type=int, default=15)
    p.add_argument(
        "--bio", action="store_true", help="emit the chromosome's attributions as a BioLang program"
    )
    p.add_argument("--out", help="with --bio: write the program here instead of printing it")
    p.set_defaults(fn=cmd_budget)

    p = sub.add_parser(
        "closure", help="the 98%%: do the attributions reproduce what each cell makes (gene level)"
    )
    p.add_argument("--chrom", required=True)
    p.add_argument("--cells", help="comma list of ENCODE cell lines with a reader and RNA-seq (default four)")
    p.add_argument("--signal", type=float, default=0.05, help="RNA-seq signal a base must reach")
    p.add_argument("--cached", action="store_true", help="show the saved result instead of recomputing")
    p.add_argument("--top", type=int, default=12)
    p.set_defaults(fn=cmd_closure)

    p = sub.add_parser(
        "decompile",
        help="one gene as a decompiled locus: every layer read for it, evidence per line, the unknown",
    )
    p.add_argument("gene")
    p.add_argument("--chrom", required=True)
    p.add_argument("--json", action="store_true", help="the layers as JSON instead of the view")
    p.set_defaults(fn=cmd_decompile)

    p = sub.add_parser(
        "motifs", help="JASPAR motifs over promoters and elements: requires lists, enrichment, operators"
    )
    p.add_argument("--chrom", help="one chromosome (sequence and GENCODE models local)")
    p.add_argument("--gene", help="with --chrom: one gene's requires list")
    p.add_argument("--library", help="one BioLib library from the genome summary")
    p.add_argument("--cached", action="store_true", help="show the saved result instead of rescanning")
    p.add_argument("--top", type=int, default=12)
    p.set_defaults(fn=cmd_motifs)

    p = sub.add_parser(
        "duplications", help="segmental duplications over the UNKNOWN blocks: curated copy-and-paste (UCSC)"
    )
    p.add_argument("--chrom", required=True)
    p.add_argument("--cached", action="store_true", help="show the saved result instead of refetching")
    p.add_argument("--top", type=int, default=10)
    p.set_defaults(fn=cmd_duplications)

    p = sub.add_parser(
        "across", help="one locus across species: alignment, conserved grammar, both constraint axes"
    )
    p.add_argument("--locus", default="ZRS", help="ZRS (the SHH limb enhancer) or HERC2_OCA2 (eye colour)")
    p.add_argument("--cached", action="store_true", help="show the saved result instead of refetching")
    p.add_argument("--top", type=int, default=15)
    p.set_defaults(fn=cmd_across)

    p = sub.add_parser(
        "profiling",
        help="phylogenetic profiling between libraries: gained and lost together across the species",
    )
    p.add_argument("--library", help="one BioLib library: its gain curve and most similar libraries")
    p.add_argument("--cached", action="store_true", help="show the saved result instead of recomputing")
    p.add_argument("--top", type=int, default=12)
    p.set_defaults(fn=cmd_profiling)

    p = sub.add_parser(
        "origin",
        help="origin per gene and age per library: the deepest clade with an orthologue (Ensembl Compara)",
    )
    p.add_argument("--gene", help="one gene symbol")
    p.add_argument("--library", help="one BioLib library id, e.g. core.translation")
    p.add_argument("--top", type=int, default=15)
    p.set_defaults(fn=cmd_origin)

    p = sub.add_parser(
        "variation",
        help="constraint on two axes, 241 mammals and 76,156 people, per UNKNOWN block and element",
    )
    p.add_argument("--chrom", help="one chromosome with a budget result")
    p.add_argument(
        "--vista-genome", action="store_true", help="VISTA positives against negatives, all chromosomes"
    )
    p.add_argument("--cached", action="store_true", help="show the saved result instead of recomputing")
    p.add_argument("--top", type=int, default=12)
    p.set_defaults(fn=cmd_variation)

    p = sub.add_parser(
        "syntax", help="syntax against values at a gene: constrained bases, where people differ, the overlap"
    )
    p.add_argument("--gene", required=True)
    p.add_argument("--chrom", required=True)
    p.add_argument("--name", help="comma list of imported individuals (default: everyone imported)")
    p.add_argument("--flank", type=int, default=0, help="bases to add on both sides of the gene")
    p.add_argument("--top", type=int, default=12)
    p.add_argument(
        "--save",
        action="store_true",
        help="save syntax_<GENE>.json under data/results (open-consent GIAB people) or the person's dir",
    )
    p.set_defaults(fn=cmd_syntax)

    p = sub.add_parser(
        "organise",
        help="the 98%%: every UNKNOWN block with both constraint axes and the copy flag, copies first",
    )
    p.add_argument("--chrom", help="one chromosome; omit to distil the genome from the saved results")
    p.add_argument("--cached", action="store_true", help="show the saved result instead of recomputing")
    p.add_argument("--blocks", action="store_true", help="list every block, not only the candidates")
    p.add_argument("--top", type=int, default=15)
    p.set_defaults(fn=cmd_organise)

    p = sub.add_parser(
        "protein", help="a gene's protein: our translation, UniProt record, AlphaFold structure"
    )
    p.add_argument(
        "--compile",
        action="store_true",
        help="compile the full definition: Ensembl, UniProt, InterPro, PDB, AlphaFold, Reactome, STRING, HPA",
    )
    p.add_argument("--refresh", action="store_true", help="ignore the local knowledge cache")
    p.add_argument(
        "--bio", action="store_true", help="emit the compiled definition as a BioLang protein block"
    )
    p.add_argument("--lib", action="store_true", help="read from the packaged proteome library (offline)")
    p.add_argument("symbol")
    p.add_argument("--chrom", help="chromosome for our own translation (e.g. chr21)")
    p.add_argument("--genome", default="data/reference/chr21.fa.gz")
    p.add_argument("--no-structure", action="store_true")
    p.set_defaults(fn=cmd_protein)

    p = sub.add_parser(
        "flow", help="walk one gene upward: DNA → regulation → RNA → protein → cell → tissue → organism"
    )
    p.add_argument("symbol")
    p.add_argument("--chrom", default="chr21")
    p.add_argument("--genome", default="data/reference/chr21.fa.gz")
    p.add_argument("--gff3")
    p.set_defaults(fn=cmd_flow)

    p = sub.add_parser("pathway", help="run a Reactome pathway as a reachability graph; knock a protein out")
    p.add_argument("ids", nargs="*", help="Reactome pathway ids, e.g. R-HSA-69541")
    p.add_argument("--gene", help="use every pathway of this gene's compiled definition")
    p.add_argument("--knockout", help="gene symbol or UniProt accession to remove")
    p.add_argument("--limit", type=int, default=50)
    p.add_argument("-v", "--verbose", action="store_true")
    p.add_argument(
        "--kinetic", action="store_true", help="run a curated BioModels ODE model of the pathway instead"
    )
    p.add_argument("--model", help="a BioModels id (BIOMD…) to run; implies --kinetic")
    p.add_argument("--query", help="search BioModels with this text instead of the pathway's name")
    p.add_argument("--hours", type=float, default=100.0, help="time units to simulate (the model's own)")
    p.add_argument("--dt", type=float, default=0.01)
    p.add_argument("--top", type=int, default=12)
    p.set_defaults(fn=cmd_pathway)

    p = sub.add_parser("proteome", help="compile every protein of a chromosome; keep the coverage summary")
    p.add_argument("--chrom", default="chr21")
    p.add_argument("--gff3")
    p.add_argument("--limit", type=int, help="first N genes only (no result saved)")
    p.set_defaults(fn=cmd_proteome)

    p = sub.add_parser(
        "segments", help="segment parser v1: genes from learned signals alone, scored vs GENCODE"
    )
    p.add_argument("--chrom", default="chr21")
    p.add_argument("--genome", default="data/reference/chr21.fa.gz")
    p.add_argument("--signals", default="data/results/signals_chr21.json")
    p.add_argument("--min-relative", type=float, default=0.6)
    p.add_argument(
        "--predicted-sites",
        action="store_true",
        help="feature c: take donors and acceptors from AlphaGenome's splice-site tracks (needs the key)",
    )
    p.add_argument(
        "--rna-filter",
        type=float,
        nargs="?",
        const=0.3,
        default=None,
        metavar="FRACTION",
        help="keep only candidates whose exons AlphaGenome predicts transcribed (mean fraction, default 0.3)",
    )
    p.add_argument(
        "--rna-measured",
        metavar="CELL_TYPE",
        help="keep candidates whose exons carry ENCODE total RNA-seq signal in these cell lines (comma list)",
    )
    p.add_argument(
        "--rna-combine",
        choices=("union", "intersection"),
        default="union",
        help="with both RNA filters: keep what either says is transcribed (union) or what both do",
    )
    p.add_argument(
        "--coding",
        choices=("codon", "markov"),
        default="codon",
        help="coding model: codon usage log-odds (default) or a 3-periodic fifth-order Markov model",
    )
    p.add_argument(
        "--promoter-anchored",
        type=int,
        nargs="?",
        const=5_000,
        default=None,
        metavar="REACH",
        help="only begin a gene within REACH bp downstream of an ENCODE promoter-like element (default 5000)",
    )
    p.set_defaults(fn=cmd_segments)

    p = sub.add_parser(
        "reader", help="reader v1: which nodes and genes a cell type reads (ENCODE DNase peaks)"
    )
    p.add_argument("--cell-type", required=True, help="ENCODE biosample term, e.g. K562, HepG2, liver")
    p.add_argument("--versus", nargs="*", help="other cell types to compare with")
    p.add_argument("--chrom", default="chr21")
    p.add_argument("--genome", default="data/reference/chr21.fa.gz")
    p.set_defaults(fn=cmd_reader)

    p = sub.add_parser("individual", help="a person's genome as a local individual: import a VCF once")
    isub = p.add_subparsers(dest="action")
    q = isub.add_parser(
        "import", help="split a VCF (GRCh38, plain or .gz, file or URL) into per-chromosome files"
    )
    q.add_argument("vcf", help="a local file or an http(s) URL, streamed once")
    q.add_argument("--name", required=True, help="how this person is called in every layer")
    q.add_argument("--note", help="where the calls come from (caller, date, consent)")
    q.add_argument("--replace", action="store_true")
    isub.add_parser("list", help="the test human and every imported person")
    q = isub.add_parser("remove", help="delete an imported person's files")
    q.add_argument("name")
    q = isub.add_parser("report", help="one Markdown page from what has been computed for the person")
    q.add_argument("--name", required=True)
    q.add_argument("--out", help="write to a file instead of printing")
    q = isub.add_parser(
        "diff", help="two people, or one against the reference: protein-changing variants by gene, as a diff"
    )
    q.add_argument("--name", required=True)
    q.add_argument("--against", default="reference", help="another local individual, or `reference`")
    q.add_argument("--top", type=int, default=40, help="genes to show")
    q.add_argument("--out", help="write the Markdown to a file instead of printing")
    q = isub.add_parser(
        "protein", help="which protein each tissue makes in this person (dominant isoform + variants)"
    )
    q.add_argument("--name", required=True)
    q.add_argument("--gene", required=True)
    q.add_argument("--chrom", required=True)
    q = isub.add_parser("regions", help="the regions a person's calls are trusted in (BED, file or URL)")
    q.add_argument("--name", required=True)
    q.add_argument("bed")
    q = isub.add_parser(
        "trio", help="Mendelian consistency of a child against both parents; de novo candidates"
    )
    q.add_argument("--name", required=True, help="the child")
    q.add_argument("--father", required=True)
    q.add_argument("--mother", required=True)
    q.add_argument("--chrom", nargs="*")
    q.add_argument("--top", type=int, default=25)
    q = isub.add_parser("coding", help="genome-wide coding SNVs by consequence and by gene (missense first)")
    q.add_argument("--name", required=True)
    q.add_argument("--chrom", nargs="*")
    q.add_argument("--top", type=int, default=25)
    q.add_argument(
        "--predict",
        action="store_true",
        help="AlphaMissense pathogenicity per missense variant (643 MB streamed once, kept under the person)",
    )
    q = isub.add_parser(
        "knockouts", help="genome-wide truncating SNVs (nonsense, start/stop lost), homozygous first"
    )
    q.add_argument("--name", required=True)
    q.add_argument("--chrom", nargs="*")
    q.add_argument("--top", type=int, default=30)
    q = isub.add_parser("screen", help="ClinVar carrier screen: pathogenic alleles the person carries")
    q.add_argument("--name", required=True)
    q.add_argument("--chrom", nargs="*", help="default: every chromosome the person has")
    q.add_argument("--top", type=int, default=30)
    q = isub.add_parser("check", help="apply the person's variants to the local reference; assembly verdict")
    q.add_argument("--name", required=True)
    q.add_argument("--chrom", nargs="*", help="default: every chromosome with rows and a local sequence")
    q = isub.add_parser("predict", help="feature d: this person's regulatory variants on one gene, predicted")
    q.add_argument("--name", required=True)
    q.add_argument("--gene", required=True)
    q.add_argument("--chrom", required=True)
    q.add_argument("--max", type=int, default=60, help="variants to score (about 8 s each, cached)")
    q.add_argument("--top", type=int, default=20)
    q.add_argument("--json", action="store_true")
    q = isub.add_parser("genes", help="one chromosome gene by gene: variants and coding consequences")
    q.add_argument("--name", default="HG002")
    q.add_argument("--chrom", default="chr21")
    q.add_argument("--top", type=int, default=25)
    q.add_argument("--json", action="store_true")
    p.set_defaults(fn=cmd_individual, action="list")

    p = sub.add_parser(
        "mouse", help="the second mammal: a mouse chromosome through the same code, nodes compared"
    )
    p.add_argument("--chrom", default="chr19")
    p.add_argument("--top", type=int, default=15)
    p.set_defaults(fn=cmd_mouse)

    p = sub.add_parser("report", help="a gene dossier: every layer about one gene in one document")
    p.add_argument("symbol")
    p.add_argument("--chrom", default="chr21")
    p.add_argument("--genome", default="data/reference/chr21.fa.gz")
    p.add_argument("--out", help="write to a file (.md or .json)")
    p.add_argument("--json", action="store_true")
    p.set_defaults(fn=cmd_report)

    p = sub.add_parser(
        "verify", help="translate every compiled gene of a chromosome and compare with UniProt"
    )
    p.add_argument("--chrom", default="chr21")
    p.add_argument("--genome", default="data/reference/chr21.fa.gz")
    p.add_argument("--top", type=int, default=15)
    p.set_defaults(fn=cmd_verify)

    p = sub.add_parser(
        "lookup",
        help="one variant through every layer: VEP, ClinVar, gnomAD, COSMIC, local trace, protein, pathways",
    )
    p.add_argument("variant", nargs="+", help="chr21:25897620 C>T (1-based, plus strand)")
    p.add_argument("--no-pathways", action="store_true")
    p.add_argument("--json", action="store_true")
    p.set_defaults(fn=cmd_lookup)

    p = sub.add_parser(
        "predict",
        help="optional AlphaGenome feature: predicted expression change per tissue for a variant",
    )
    p.add_argument("variant", nargs="*", help="chr21:25897620 C>T (1-based, plus strand)")
    p.add_argument("--element", help="feature b: delete a region (chr21:START-END) and read which gene moves")
    p.add_argument("--status", action="store_true", help="whether the feature is enabled and what it brings")
    p.add_argument("--top", type=int, default=20)
    p.add_argument("--min", type=float, default=0.05, help="smallest |log2 fold change| reported")
    p.add_argument("--json", action="store_true")
    p.set_defaults(fn=cmd_predict)

    p = sub.add_parser("ptm", help="post-translational state: modifiable sites and who writes them")
    p.add_argument("protein", nargs="?", help="one protein's sites; omit for the genome-wide summary")
    p.add_argument("--writer", help="substrates of a kinase, acetyltransferase, … (UniProt 'by')")
    p.add_argument("--top", type=int, default=25)
    p.add_argument("--rebuild", action="store_true", help="rebuild the index from the cached definitions")
    p.set_defaults(fn=cmd_ptm)

    p = sub.add_parser("graph", help="the local protein knowledge graph built from compiled definitions")
    p.add_argument("gene", nargs="?", help="show one protein's neighbourhood instead of the summary")
    p.add_argument("--min-score", type=float, default=0.7)
    p.add_argument("--max-nodes", type=int, default=40)
    p.add_argument("--top", type=int, default=25)
    p.add_argument("--save", help="save the summary as data/results/graph_<name>.json")
    p.set_defaults(fn=cmd_graph)

    p = sub.add_parser(
        "rna", help="the RNA layer of a gene: transcripts and isoforms, GTEx expression per tissue"
    )
    p.add_argument("symbol")
    p.add_argument("--chrom", default="chr21")
    p.add_argument("--genome", default="data/reference/chr21.fa.gz")
    p.add_argument("--top", type=int, default=15)
    p.add_argument("--no-expression", action="store_true")
    p.add_argument("--isoforms", action="store_true", help="per-transcript expression per tissue (GTEx)")
    p.set_defaults(fn=cmd_rna)

    p = sub.add_parser("repeats", help="curated repeat annotation (RepeatMasker via UCSC) for a chromosome")
    p.add_argument("--chrom", default="chr21")
    p.add_argument("--genome", default="data/reference/chr21.fa.gz")
    p.add_argument("--refresh", action="store_true")
    p.set_defaults(fn=cmd_repeats)

    p = sub.add_parser(
        "regulation", help="the regulatory input of a gene: promoter, enhancers in its node, insulators"
    )
    p.add_argument("--gene", required=True)
    p.add_argument("--chrom", default="chr21")
    p.add_argument("--genome", default="data/reference/chr21.fa.gz")
    p.add_argument("--gff3")
    p.add_argument("--top", type=int, default=15)
    p.set_defaults(fn=cmd_regulation)

    p = sub.add_parser("domains", help="nodes above genes: domains inferred from CTCF boundaries")
    p.add_argument("--chrom", default="chr21")
    p.add_argument("--genome", default="data/reference/chr21.fa.gz")
    p.add_argument("--gff3")
    p.add_argument("--top", type=int, default=12)
    p.add_argument(
        "--contact-map",
        action="store_true",
        help="hold the CTCF-only boundaries against insulation minima of AlphaGenome's predicted contact map",
    )
    p.add_argument(
        "--hic",
        metavar="BIOSOURCE",
        help="hold the CTCF-only boundaries against 4DN's measured boundary calls (needs a 4DN key)",
    )
    p.set_defaults(fn=cmd_domains)

    p = sub.add_parser("libs", help="list the biological libraries found in the genome")
    p.add_argument(
        "--proteome", action="store_true", help="the packaged proteome library: counts and evidence"
    )
    p.add_argument("--layer", choices=LAYERS)
    p.add_argument("-v", "--verbose", action="store_true")
    p.add_argument("--data", action="store_true", help="count members from GO/Reactome data")
    p.add_argument("--verify", action="store_true", help="check catalogue genes against the data")
    p.add_argument("--gene", help="which libraries a gene belongs to, by data")
    p.set_defaults(fn=cmd_libs)

    p = sub.add_parser(
        "therapeutic",
        help="therapeutic target and mechanism reasoning for one tumour (research hypotheses only)",
    )
    p.add_argument("--tumour", required=True, help="tumour VCF")
    p.add_argument("--normal", help="matched normal VCF; without it germline is estimated")
    p.add_argument("--rna", help="tumour RNA table: gene<TAB>value per line")
    p.add_argument("--hla", action="append", help="patient class I alleles, e.g. HLA-A*02:01,HLA-B*07:02")
    p.add_argument(
        "--hla-predictor",
        default="",
        help="peptide/HLA binding predictor: mhcflurry, netmhcpan, none (default: first installed)",
    )
    p.add_argument(
        "--install-hla-models",
        action="store_true",
        help="download MHCflurry's class I models (its own downloader fails on Python 3.13+)",
    )
    p.add_argument("--chrom", action="append", help="restrict to these chromosomes")
    p.add_argument("--purity", type=float, help="tumour purity 0-1, for clonality")
    p.add_argument("--cnv", help="copy-number table: gene<TAB>copies per line")
    p.add_argument(
        "--cohort",
        help="cBioPortal study id of the same cancer type, e.g. brca_tcga_pan_can_atlas_2018; "
        "gives how often each gene is raised above normal tissue in that cancer",
    )
    p.add_argument("--top", type=int, default=12, help="genes to analyse, by variant rank")
    p.add_argument("--detail", type=int, default=5, help="candidates to expand in the report")
    p.add_argument("--report", action="store_true", help="print the full text report")
    p.add_argument("--spec", help="print the target specification for one gene")
    p.add_argument("--out", help="directory for the dataset layers and the report")
    p.add_argument("--offline", action="store_true", help="local caches only, no network")
    p.add_argument("--no-indirect", action="store_true", help="skip pathway-induced candidates")
    p.set_defaults(fn=cmd_therapeutic)

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
    from genomeos.lib.biolang import register as _register_imports

    _register_imports()  # `import protein:SYMBOL` in BioLang programs reads the packaged proteome
    args = build_parser().parse_args(argv)
    # a command that defaults to chromosome 21's sequence follows --chrom when that chromosome is local
    chrom = getattr(args, "chrom", None)
    if isinstance(chrom, str) and getattr(args, "genome", None) == "data/reference/chr21.fa.gz":
        local = Path(f"data/reference/{chrom}.fa.gz")
        if chrom != "chr21" and local.exists():
            args.genome = str(local)
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
