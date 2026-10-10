# SPDX-License-Identifier: AGPL-3.0-or-later
"""M5: which public measurement could score milestone 1.3's clause 2, counted from coordinates alone.

    uv run python scripts/measured_unknown_inventory.py [--no-save] [--offline]

The clause asks for the 882 real-unknown blocks attributed to a gene and a tissue. The project's own
assays reach 3 of them with an element CRISPRi ever tested against a coding gene (`a7f207f`), and the
calibrated sizing (`cab70d9`) says no design reaches 80% probability of detection below 75 compared
blocks. Wet-lab synthesis is a permanent non-goal (docs/DECISIONS.md D8), so this asks the only open
question left: does anything already published reach enough of those blocks?

For every candidate source it records where the data is (URL, accession, bytes, sha256 and the
publisher's md5 where there is one), the licence, the build and the lift, and then counts, through
`genomeos.attribution.measured_unknown.census`, how many of the 882 blocks (and of the 531 carrying a
scored element) a tested interval reaches. **It reads coordinates, the gene an interval was measured
against and the cell it was measured in. It reads no effect size, activity, p-value, call or PIP**:
the parsers select columns by name, and the two range-read tracks are reduced to positions before any
count. Raw files go to data/knowledge/measured_unknown (git-ignored); the result to data/results.

No model request, no key, no money.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import Any

from genomeos import manifest as mf
from genomeos.attribution import closure_cell, crispri, human_panel, mpra, vista
from genomeos.attribution import measured_unknown as mu
from genomeos.results import save_result

ROOT = Path(__file__).resolve().parent.parent
ENTRY = "scripts/measured_unknown_inventory.py"
OWN_CODE = ("genomeos/attribution/measured_unknown.py", ENTRY, "tests/test_measured_unknown.py")
NAME = "measured_unknown_inventory"
ENCODE = "https://www.encodeproject.org"
ENCODE_LICENCE = "ENCODE data use policy (no restrictions)"
CHAIN_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120861/suppl/"
GASPERINI_PAIRS = GEO + "GSE120861_gene_gRNAgroup_pair_table.at_scale.txt.gz"
DBSNP_COMMON = f"{human_panel.GBDB}/snp/dbSnp155Common.bb"
UKB_FINEMAP = "https://www.dropbox.com/s/cdsdgwxkxkcq8cn/UKBB_94traits_release1.1.tar.gz?dl=1"
MPRABASE = "https://mprabase.ucsf.edu/app/mprabase"
BENCHMARK = (
    (
        "EPCrisprBenchmark_combined_data.training_K562.GRCh38.tsv.gz",
        "epcrispr_training_k562",
        "training, K562",
    ),
    (
        "EPCrisprBenchmark_combined_data.heldout_5_cell_types.GRCh38.tsv.gz",
        "epcrispr_heldout_5",
        "held-out, 5 cell types",
    ),
)
USER_AGENT = "GenomeOS/0.9 (M5 measured-unknown inventory; coordinates only)"
#: the guide-sequence check that decides a build when a file does not state one
ASSEMBLY_SAMPLE = 200
ASSEMBLY_AGREE = 0.9


def say(*a: Any) -> None:
    print(*a, flush=True)


# ---- fetching, with every byte digested -----------------------------------------------------------------


def http_json(url: str) -> Any:
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=180) as r:
        return json.load(r)


def fetch(url: str, name: str, md5: str | None = None, offline: bool = False) -> dict[str, Any]:
    """`url` cached at data/knowledge/measured_unknown/`name`, with its sha256, md5 and size."""
    path = mu.STORE / name
    if not path.exists():
        if offline:
            raise FileNotFoundError(path)
        mu.STORE.mkdir(parents=True, exist_ok=True)
        req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        tmp = path.with_suffix(path.suffix + ".part")
        with urllib.request.urlopen(req, timeout=600) as r, open(tmp, "wb") as out:
            while chunk := r.read(1 << 20):
                out.write(chunk)
        tmp.rename(path)
    h, m = hashlib.sha256(), hashlib.md5()  # noqa: S324 - md5 only to compare with the publisher's
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
            m.update(chunk)
    return {
        "url": url,
        "path": str(path),
        "bytes": path.stat().st_size,
        "sha256": h.hexdigest(),
        "md5": m.hexdigest(),
        "md5_published": md5,
        "md5_matches": None if md5 is None else m.hexdigest() == md5,
    }


def lines(path: str | Path) -> Iterator[str]:
    raw = Path(path).read_bytes()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    yield from io.StringIO(raw.decode("utf-8", "replace"))


def encode_files(**params: str) -> tuple[str, list[dict[str, Any]]]:
    """One ENCODE portal search for released files; the query URL is returned to be recorded."""
    q = {"type": "File", "status": "released", "format": "json", "limit": "all"}
    q.update(params)
    fields = ("accession", "dataset", "assembly", "file_format", "file_size", "md5sum", "href", "output_type")
    url = f"{ENCODE}/search/?" + urllib.parse.urlencode(q) + "".join(f"&field={f}" for f in fields)
    return url, http_json(url).get("@graph", [])


def encode_experiments(assay: str) -> tuple[str, dict[str, dict[str, Any]]]:
    """Released functional characterization experiments of one assay: biosample, gene, references."""
    fields = (
        "accession",
        "assay_title",
        "biosample_ontology.term_name",
        "examined_loci.gene.symbol",
        "elements_references.accession",
        "lab.title",
    )
    q = {
        "type": "FunctionalCharacterizationExperiment",
        "status": "released",
        "assay_title": assay,
        "format": "json",
        "limit": "all",
    }
    url = f"{ENCODE}/search/?" + urllib.parse.urlencode(q) + "".join(f"&field={f}" for f in fields)
    return url, {x["@id"]: x for x in http_json(url).get("@graph", [])}


def encode_fetch(f: dict[str, Any], offline: bool) -> dict[str, Any]:
    name = f["href"].rsplit("/", 1)[-1]
    got = fetch(ENCODE + f["href"], name, f.get("md5sum"), offline)
    got["accession"] = f["accession"]
    got["assembly"] = f.get("assembly")
    return got


def biosample(x: dict[str, Any]) -> str:
    return (x.get("biosample_ontology") or {}).get("term_name") or ""


# ---- the build of files that do not state one -------------------------------------------------------------


def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def sequence_agreement(guides: list[tuple[mu.Tested, str]]) -> dict[str, Any]:
    """The share of guide spacers found on GRCh38 at their stated coordinates (either strand, with the
    spacer's own length of slack for 0/1-based and PAM conventions, and with or without a leading
    added G). High means GRCh38."""
    from genomeos.genome.index import IndexedGenome, has_reference, reference_fasta
    from genomeos.genome.sequence import Locus, Strand

    by_chrom: dict[str, list[tuple[mu.Tested, str]]] = {}
    for t, spacer in guides:
        if len(spacer) >= 17 and set(spacer) <= set("ACGT"):
            by_chrom.setdefault(t.chrom, []).append((t, spacer))
    if not by_chrom:
        return {"checked": 0, "agree": 0, "share": None}
    chrom = max(by_chrom, key=lambda c: len(by_chrom[c]))
    if not has_reference(chrom):
        return {"checked": 0, "agree": 0, "share": None, "chrom": chrom}
    sample = by_chrom[chrom][:ASSEMBLY_SAMPLE]
    agree = 0
    with IndexedGenome(reference_fasta(chrom)) as g:
        for t, spacer in sample:
            lo = max(0, t.start - 3 - len(spacer))
            ref = str(g.fetch(Locus(chrom, lo, t.end + 3 + len(spacer), Strand.PLUS))).upper()
            # a 21-nt guide usually carries a 5' G added for the U6 promoter, absent from the genome
            if any(x in ref or revcomp(x) in ref for x in (spacer, spacer[1:])):
                agree += 1
    return {"chrom": chrom, "checked": len(sample), "agree": agree, "share": round(agree / len(sample), 4)}


def decide_build(
    chain: mu.Chain, guides: list[tuple[mu.Tested, str]]
) -> tuple[list[mu.Tested], dict[str, Any]]:
    """Guides whose file states no build: kept as GRCh38 if their spacers sit there, lifted from hg19 if
    they sit there after the lift, and dropped (recorded) otherwise."""
    as38 = sequence_agreement(guides)
    if (as38["share"] or 0) >= ASSEMBLY_AGREE:
        return [t for t, _ in guides], {"grch38": as38, "build": "GRCh38"}
    lifted = []
    for t, s in guides:
        got = chain.lift(t.chrom, t.start, t.end)
        if got:
            lifted.append((mu.Tested(got[0], got[1], got[2], t.gene, t.cell), s))
    back = sequence_agreement(lifted)
    build = "hg19" if (back["share"] or 0) >= ASSEMBLY_AGREE else "unresolved"
    check = {
        "grch38": as38,
        "after_hg19_lift": back,
        "build": build,
        "lifted": len(lifted),
        "rows": len(guides),
    }
    return ([t for t, _ in lifted] if build == "hg19" else []), check


# ---- the sources ----------------------------------------------------------------------------------------


def source(
    key: str,
    label: str,
    endpoint: str,
    rows: list[mu.Tested] | None,
    *,
    held: bool,
    cells: Iterable[str],
    licence: str,
    files: list[dict[str, Any]],
    assembly: str,
    query: str | None = None,
    lifted: dict[str, Any] | None = None,
    gene_resolved: bool = False,
    note: str = "",
    not_computed: str | None = None,
) -> dict[str, Any]:
    return {
        "key": key,
        "label": label,
        "endpoint": endpoint,
        "endpoint_meaning": mu.ENDPOINTS[endpoint],
        "held_before_this_lane": held,
        "cell_types": sorted(set(cells) - {""}),
        "licence": licence,
        "assembly": assembly,
        "lift": lifted,
        "query": query,
        "files": files,
        "bytes": sum(f.get("bytes") or 0 for f in files),
        "gene_resolved_per_interval": gene_resolved,
        "note": note,
        "not_computed": not_computed,
        "_rows": rows,
    }


def held_benchmark(name: str, key: str, label: str) -> dict[str, Any]:
    pairs = crispri.load(name)
    rows = [mu.Tested(p.chrom, p.start, p.end, p.gene, p.cell) for p in pairs]
    return source(
        key,
        f"EPCrisprBenchmark (ENCODE-rE2G) {label}, as held",
        "gene",
        rows,
        held=True,
        cells=[p.cell for p in pairs],
        licence="MIT (EngreitzLab/CRISPR_comparison repository licence); Gschwind et al. 2025",
        files=[
            {"url": crispri.BASE_URL + name, **mf.input_entry(crispri.KNOWLEDGE / name, partition="held")}
        ],
        assembly="GRCh38",
        gene_resolved=True,
        note="the ENCODE-rE2G benchmark set: element-gene pairs after the benchmark's own filters",
    )


def flowfish(chain: mu.Chain, offline: bool) -> list[dict[str, Any]]:
    """ENCODE CRISPRi-FlowFISH: the series' element quantifications (element, measured gene), and the
    libraries' guide targets, which bound what every experiment could have measured (gene unresolved)."""
    q_exp, exps = encode_experiments("Flow-FISH CRISPR screen")
    q, files = encode_files(assay_title="Flow-FISH CRISPR screen", output_type="element quantifications")
    cols = ("chrPerturbationTarget", "startPerturbationTarget", "endPerturbationTarget")
    header: list[str] | None = None
    got_files, hg19_rows, hg38_rows = [], [], []
    for f in sorted(files, key=lambda f: f["accession"]):
        if f.get("file_format") == "tsv" and f.get("assembly") == "hg19":
            got = encode_fetch(f, offline)
            text = list(lines(got["path"]))
            header = text[0].rstrip("\n").split("\t")
            hg19_rows += list(mu.table_rows(text, *cols, gene="measuredGeneSymbol"))
            got_files.append(got)
    for f in sorted(files, key=lambda f: f["accession"]):
        if f.get("file_format") == "bed" and f.get("assembly") == "GRCh38" and header:
            got = encode_fetch(f, offline)
            # the GRCh38 BED is the same table without its header line: the tsv's header names it
            body = [ln for ln in lines(got["path"]) if ln.strip()]
            hg38_rows += list(
                mu.table_rows(["\t".join(header) + "\n", *body], *cols, gene="measuredGeneSymbol")
            )
            got_files.append(got)
    lifted, lost = mu.lift_all(chain, hg19_rows)
    cells = [biosample(x) for x in exps.values()]
    out = [
        source(
            "encode_flowfish_elements",
            "ENCODE CRISPRi-FlowFISH, series element quantifications (Engreitz, Sabeti labs)",
            "gene",
            lifted + hg38_rows,
            held=False,
            cells=cells,
            licence=ENCODE_LICENCE,
            files=got_files,
            assembly="hg19 tsv lifted to GRCh38; GRCh38 bed as published",
            query=q,
            lifted={
                "rows_hg19": len(hg19_rows),
                "lifted": len(lifted),
                "not_lifted": lost,
                "rows_grch38": len(hg38_rows),
            },
            gene_resolved=True,
            note=f"{len(exps)} released experiments ({q_exp}); element-gene pairs as tested, before any "
            "benchmark filter",
        )
    ]
    refs = Counter(r["accession"] for x in exps.values() for r in x.get("elements_references", []))
    _, rfiles = encode_files(output_type="elements reference")
    lib_rows: list[mu.Tested] = []
    lib_files: list[dict[str, Any]] = []
    checks: dict[str, Any] = {}
    for f in sorted(rfiles, key=lambda f: f["accession"]):
        ref = f["dataset"].strip("/").split("/")[-1]
        if ref not in refs or f.get("file_format") != "tsv":
            continue
        got = encode_fetch(f, offline)
        lib_files.append(got)
        text = list(lines(got["path"]))
        head = text[0].rstrip("\n").split("\t") if text else []
        if head[:2] == ["OligoID", "target"]:
            rows = list(mu.region_column_rows(text, "target"))
            if f.get("assembly") == "hg19":
                rows, lost = mu.lift_all(chain, rows)
                checks[f["accession"]] = {"build": "hg19 as stated", "not_lifted": lost}
            else:
                checks[f["accession"]] = {"build": "unstated; dropped"}
                rows = []
            lib_rows += rows
        elif {"chr", "start", "end", "GuideSequence"} <= set(head):
            ic, i0, i1, ig = (head.index(c) for c in ("chr", "start", "end", "GuideSequence"))
            guides = []
            for ln in text[1:]:
                v = ln.rstrip("\n").split("\t")
                try:
                    t = mu.Tested(v[ic], int(float(v[i0])), int(float(v[i1])))
                except (ValueError, IndexError):
                    continue
                if t.chrom.startswith("chr") and t.end > t.start:
                    guides.append((t, v[ig].upper()))
            rows, checks[f["accession"]] = decide_build(chain, guides)
            lib_rows += rows
        else:
            checks[f["accession"]] = {"build": "no coordinate columns; a region list, not read as guides"}
    out.append(
        source(
            "encode_flowfish_library_targets",
            "ENCODE CRISPRi-FlowFISH guide libraries (elements references), gene not resolved per guide",
            "gene",
            lib_rows,
            held=False,
            cells=cells,
            licence=ENCODE_LICENCE,
            files=lib_files,
            assembly="per file: stated hg19 lifted; unstated decided by guide sequence against GRCh38",
            lifted={"per_file": checks},
            note="an upper bound: a library target was measured against the gene of at least one "
            "experiment using the library, not necessarily the block's attributed gene",
        )
    )
    return out


def gasperini(chain: mu.Chain, offline: bool) -> dict[str, Any]:
    got = fetch(GASPERINI_PAIRS, "GSE120861_gene_gRNAgroup_pair_table.at_scale.txt.gz", None, offline)
    rows = list(
        mu.table_rows(
            lines(got["path"]),
            "gRNAgroup.chr",
            "gRNAgroup.start",
            "gRNAgroup.stop",
            gene="targetgene_short_name",
            cell="K562",
        )
    )
    lifted, lost = mu.lift_all(chain, rows)
    return source(
        "gasperini2019_at_scale",
        "Gasperini et al. 2019, CRISPRi with single-cell RNA-seq, at-scale screen (GSE120861)",
        "gene",
        lifted,
        held=False,
        cells=["K562"],
        licence="NCBI GEO public record GSE120861 (no licence stated; cite Gasperini et al. 2019, Cell)",
        files=[got],
        assembly="hg19, lifted to GRCh38",
        lifted={"rows_hg19": len(rows), "lifted": len(lifted), "not_lifted": lost},
        gene_resolved=True,
        note="every gRNA-group by gene pair the screen tested (the pair table); the differential "
        "expression table that holds the readout is not fetched",
    )


def gasperini_groups(offline: bool) -> set[str]:
    """The at-scale screen's gRNA-group identifiers without their suffix ('chr12.3641_top_two')."""
    got = fetch(GASPERINI_PAIRS, "GSE120861_gene_gRNAgroup_pair_table.at_scale.txt.gz", None, offline)
    it = lines(got["path"])
    head = next(it).rstrip("\n").split("\t")
    i = head.index("gRNAgroup")
    return {ln.split("\t")[i].split("_", 1)[0] for ln in it if ln.strip()}


def scrna_guides(chain: mu.Chain, offline: bool) -> dict[str, Any]:
    """ENCODE 'perturbation followed by scRNA-seq' guide libraries: whole transcriptome, gene unresolved."""
    q, exps = encode_experiments("perturbation followed by scRNA-seq")
    refs = {r["accession"] for x in exps.values() for r in x.get("elements_references", [])}
    _, rfiles = encode_files(output_type="elements reference", file_format="fasta")
    rows: list[mu.Tested] = []
    files, checks = [], {}
    for f in sorted(rfiles, key=lambda f: f["accession"]):
        if f["dataset"].strip("/").split("/")[-1] not in refs:
            continue
        got = encode_fetch(f, offline)
        files.append(got)
        guides = [(t, s) for t, _, s in mu.fasta_name_rows(lines(got["path"]))]
        if guides:
            kept, checks[f["accession"]] = decide_build(chain, guides)
            rows += kept
            continue
        # no coordinates in the names: say whose library it is, by its group identifiers
        ids = {ln[1:].split("_", 1)[0] for ln in lines(got["path"]) if ln.startswith(">chr")}
        checks[f["accession"]] = {
            "build": "no coordinates in the guide names",
            "group_ids": len(ids),
            "ids_in_gasperini2019_at_scale": len(ids & gasperini_groups(offline)),
            "counted_under": "gasperini2019_at_scale",
        }
    return source(
        "encode_perturb_scrnaseq_guides",
        "ENCODE CRISPRi with single-cell RNA-seq guide libraries (Shendure, Reddy labs)",
        "gene",
        rows,
        held=False,
        cells=[biosample(x) for x in exps.values()],
        licence=ENCODE_LICENCE,
        files=files,
        assembly="decided per file by guide sequence against GRCh38",
        query=q,
        lifted={"per_file": checks},
        note="whole-transcriptome readout: a guide was measured against every expressed gene, so the "
        "gene is not resolved per guide; an upper bound",
    )


def facs(offline: bool) -> dict[str, Any]:
    """ENCODE FACS CRISPR screens on human cells (CREST-seq and kin): one guide table per experiment."""
    q_exp, exps = encode_experiments("FACS CRISPR screen")
    q, files = encode_files(assay_title="FACS CRISPR screen", output_type="guide quantifications")
    per_exp: dict[str, dict[str, Any]] = {}
    for f in sorted(files, key=lambda f: f["accession"]):
        if f.get("assembly") == "GRCh38" and f["dataset"] not in per_exp:
            per_exp[f["dataset"]] = f
    rows: list[mu.Tested] = []
    got_files = []
    for ds, f in sorted(per_exp.items()):
        x = exps.get(ds, {})
        genes = sorted(
            {(locus.get("gene") or {}).get("symbol") or "" for locus in x.get("examined_loci", [])}
        )
        got = encode_fetch(f, offline)
        got_files.append(got)
        for t in mu.table_or_bed_rows(
            list(lines(got["path"])),
            "chrPerturbationTarget",
            "startPerturbationTarget",
            "endPerturbationTarget",
        ):
            for g in genes or [""]:
                rows.append(mu.Tested(t.chrom, t.start, t.end, g or None, biosample(x)))
    return source(
        "encode_facs_crispr_human",
        "ENCODE FACS CRISPR screens, human cells (CREST-seq tiling; one sorted reporter gene each)",
        "gene",
        rows,
        held=False,
        cells=[biosample(exps.get(ds, {})) for ds in per_exp],
        licence=ENCODE_LICENCE,
        files=got_files,
        assembly="GRCh38 as published",
        query=q,
        gene_resolved=True,
        note=f"one GRCh38 guide table per experiment ({len(per_exp)}; {q_exp}); the gene is the "
        "experiment's examined locus",
    )


def cell_level(offline: bool) -> dict[str, Any]:
    """Cell-growth screens: the held whole-genome K562 guide map and the published element tables."""
    rows = [
        mu.Tested(c, s, e, None, "K562")
        for c, v in closure_cell.guides_by_chrom(closure_cell.GUIDE_BED).items()
        for s, e, _ in v
    ]
    files = [mf.input_entry(closure_cell.GUIDE_BED, partition="held (lane-closurecell)")]
    queries = []
    for assay in ("proliferation CRISPR screen", "CRISPR screen"):
        q, fl = encode_files(assay_title=assay, output_type="element quantifications", assembly="GRCh38")
        queries.append(q)
        for f in sorted(fl, key=lambda f: f["accession"]):
            got = encode_fetch(f, offline)
            files.append(got)
            text = list(lines(got["path"]))
            try:
                rows += mu.table_or_bed_rows(
                    text, "chrPerturbationTarget", "startPerturbationTarget", "endPerturbationTarget"
                )
            except KeyError:
                rows += list(mu.bed_rows(text))
    return source(
        "encode_cell_growth_screens",
        "ENCODE CRISPR growth and survival screens (whole-genome K562 guide map; element tables)",
        "cell",
        rows,
        held=True,
        cells=["K562", "WTC11", "HCT116", "others listed on the portal"],
        licence=ENCODE_LICENCE,
        files=files,
        assembly="GRCh38 as published",
        query=" ; ".join(queries),
        note="the readout is the cell's growth, not a gene's expression (lane-closurecell, M3)",
    )


def mpra_sources(chain: mu.Chain, offline: bool) -> list[dict[str, Any]]:
    held_rows = [mu.Tested(e.chrom, e.start, e.end) for e in mpra.load()]
    out = [
        source(
            "encode_lentimpra_joint_library",
            "ENCODE4 lentiMPRA joint library ENCSR106SZM (K562, HepG2, WTC11), as held",
            "activity",
            held_rows,
            held=True,
            cells=["K562", "HepG2", "WTC11"],
            licence=ENCODE_LICENCE,
            files=[
                {
                    "url": mpra.ENCODE_FILE.format(acc=a),
                    **mf.input_entry(mpra.KNOWLEDGE / f"{a}.bed.gz", partition="held"),
                }
                for a in mpra.FILES.values()
            ],
            assembly="GRCh38 as published",
            note="Agarwal et al. 2025; the project's measured layer since 2026-09-16",
        )
    ]
    q_exp, exps = encode_experiments("MPRA")
    refs: dict[str, set[str]] = {}
    for x in exps.values():
        for r in x.get("elements_references", []):
            refs.setdefault(r["accession"], set()).add(biosample(x))
    _, rfiles = encode_files(output_type="elements reference")
    rows: list[mu.Tested] = []
    files, cells, skipped = [], set(), []
    for f in sorted(rfiles, key=lambda f: f["accession"]):
        ref = f["dataset"].strip("/").split("/")[-1]
        if ref not in refs or ref == mpra.LIBRARY:
            continue
        fmt, asm = f.get("file_format"), f.get("assembly")
        if fmt not in ("bed", "tsv") or asm not in ("GRCh38", "hg19"):
            skipped.append({"accession": f["accession"], "format": fmt, "assembly": asm, "reference": ref})
            continue
        got = encode_fetch(f, offline)
        files.append(got)
        cells |= refs[ref]
        text = list(lines(got["path"]))
        if fmt == "bed":
            got_rows = list(mu.bed_rows(text))
        else:
            at = next(
                (
                    i
                    for i, ln in enumerate(text)
                    if ln.split("\t")[:5] == ["name", "type", "chr", "start", "stop"]
                ),
                None,
            )
            got_rows = [] if at is None else list(mu.table_rows(text[at:], "chr", "start", "stop"))
            if at is None:
                skipped.append({"accession": f["accession"], "why": "no name/type/chr/start/stop header"})
        if asm == "hg19":
            got_rows, _ = mu.lift_all(chain, got_rows)
        rows += got_rows
    out.append(
        source(
            "encode_mpra_element_references",
            "ENCODE MPRA and lentiMPRA libraries other than the joint library (elements references)",
            "activity",
            rows,
            held=False,
            cells=cells,
            licence=ENCODE_LICENCE,
            files=files,
            assembly="GRCh38 as published; hg19 lifted",
            query=q_exp,
            note="includes the large cell-type lentiMPRA libraries of Agarwal et al. 2025; the variant "
            "libraries published only as FASTA are listed under not_computed",
            not_computed=f"{len(skipped)} reference files carry no coordinate table (FASTA or no build): "
            + json.dumps(skipped[:8]),
        )
    )
    return out


def vista_source() -> dict[str, Any]:
    rows = [mu.Tested(e.chrom, e.start, e.end) for e in vista.load_loci()]
    return source(
        "vista",
        "VISTA Enhancer Browser, transgenic mouse e11.5 (LBNL)",
        "activity",
        rows,
        held=True,
        cells=["mouse embryo e11.5"],
        licence="VISTA Enhancer Browser, public with citation (Visel et al. 2007)",
        files=[
            {
                "url": vista.VISTA_LOCUS_URL,
                **mf.input_entry(vista.KNOWLEDGE / "locus.tsv.gz", partition="held"),
            }
        ],
        assembly="GRCh38 as held",
        note="positive and negative tests alike",
    )


def track_positions(url: str, rows: list[dict[str, Any]]) -> tuple[list[mu.Tested], dict[str, Any]]:
    """A remote bigBed reduced to item positions over the blocks; no other field is kept."""
    bb = human_panel.BigBed(url)
    r0, b0 = bb.src.requests, bb.src.bytes_fetched
    out: list[mu.Tested] = []
    by: dict[str, list[tuple[int, int]]] = {}
    for r in rows:
        by.setdefault(r["chrom"], []).append((r["start"], r["end"]))
    for chrom, ivs in sorted(by.items()):
        for item in bb.query(chrom, human_panel.merge_intervals(ivs)):
            s, e = int(item["chromStart"]), int(item["chromEnd"])
            out.append(mu.Tested(chrom, s, max(e, s + 1)))
    return out, {
        "url": url,
        "requests": bb.src.requests - r0,
        "bytes": bb.src.bytes_fetched - b0,
        "read_by": "HTTP range requests over the 882 blocks only; never downloaded whole",
    }


def mpravardb(rows: list[dict[str, Any]]) -> dict[str, Any]:
    pos, cost = track_positions(human_panel.BIGBEDS["mpravardb"], rows)
    return source(
        "mpravardb",
        "MPRAVarDB allele pairs tested in MPRAs (UCSC track mpravardb)",
        "activity",
        pos,
        held=False,
        cells=["many, per row"],
        licence="MPRAVarDB (Jin et al. 2024) as redistributed by UCSC; terms not checked beyond that",
        files=[cost],
        assembly="GRCh38 (UCSC hg38 track)",
        note="variant-level reporter tests; every column but the position is dropped at read",
    )


def gtex_tested_universe(rows: list[dict[str, Any]]) -> dict[str, Any]:
    pos, cost = track_positions(DBSNP_COMMON, rows)
    snv = [t for t in pos if t.end - t.start == 1]
    return source(
        "gtex_v8_finemapped_eqtl",
        "GTEx v8 cis-eQTL fine-mapping (DAP-G; UCSC gtexDapg), tested universe by proxy",
        "variant_gene",
        snv,
        held=True,
        cells=["49 GTEx tissues"],
        licence="GTEx portal open access (summary statistics); dbSNP public domain",
        files=[cost],
        assembly="GRCh38",
        note="coverage counted on a PROXY of the tested universe: dbSNP 155 common single-nucleotide "
        "variants (minor allele 1% or more in a major population). GTEx tested variants with minor "
        "allele 1% or more among its 838 donors within 1 Mb of a gene. The fine-mapped track itself "
        "is the readout and is not read here",
    )


def not_fetched() -> list[dict[str, Any]]:
    return [
        source(
            "encode_starr_seq",
            "ENCODE STARR-seq (whole-genome STARR-seq in K562 and A549; ATAC-STARR; others)",
            "activity",
            None,
            held=False,
            cells=["K562", "A549", "HepG2", "MCF-7", "SH-SY5Y", "T-helper 17", "neurons", "pancreas"],
            licence=ENCODE_LICENCE,
            files=[],
            assembly="GRCh38",
            not_computed="the tested universe of a whole-genome library is the genome at a depth only the "
            "input-library signal states (bigWig, 80 files, 42.9 GB); the element files are called "
            "regions, that is readouts, and are not read here. At any depth no gene is named",
        ),
        source(
            "mprabase",
            "MPRAbase (Zhao et al. 2025, Genome Research): 130 experiments, 17.7 M elements, 35 cell types",
            "activity",
            None,
            held=False,
            cells=[],
            licence="not stated where checked",
            files=[{"url": MPRABASE, "bytes": 0}],
            assembly="per experiment",
            not_computed="an interactive app with no static bulk download found; its TLS chain did not "
            "verify from this host. Its ENCODE-hosted experiments are counted above; no gene is named",
        ),
        source(
            "ukb_finemapped",
            "UK Biobank fine-mapping, 94 traits (Kanai et al. 2021, release 1.1)",
            "trait",
            None,
            held=False,
            cells=[],
            licence="not stated on the release page (finucanelab.org/data)",
            files=[{"url": UKB_FINEMAP, "bytes": 889_779_304}],
            assembly="GRCh37 and GRCh38 per the release",
            not_computed="not fetched (890 MB): the endpoint names a trait, not a gene or a tissue, and the "
            "tested universe is the common-variant universe counted under the GTEx row",
        ),
    ]


# ---- the run ----------------------------------------------------------------------------------------------


def coding_symbols() -> set[str]:
    from genomeos.genome import Annotation
    from genomeos.genome.annotation import default_gencode

    out: set[str] = set()
    for c in mu.CHROMS:
        path = default_gencode({c})
        if path is None:
            continue
        out |= {g.symbol for g in Annotation.from_gff3(path, {c}).protein_coding() if g.symbol}
    return out


def gene_union(rows: list[dict[str, Any]], sources: list[dict[str, Any]], coding: set[str]) -> dict[str, Any]:
    """Every gene-level perturbation source pooled: the clause's best case from what is published."""
    gene = [s for s in sources if s["endpoint"] == "gene" and s["_rows"] is not None]
    resolved = [t for s in gene if s["gene_resolved_per_interval"] for t in s["_rows"]]
    upper = [t for s in gene for t in s["_rows"]]
    a = mu.census(rows, resolved, coding)
    b = mu.census(rows, upper, None)
    return {
        "sources": [s["key"] for s in gene],
        "gene_resolved_sources": a,
        "with_unresolved_libraries_as_upper_bound": {
            k: b[k]
            for k in (
                "blocks_touched",
                "blocks_touched_carrying",
                "blocks_element_matched",
                "tested_inside",
                "blocks_with_tested_inside",
                "touched_blocks",
            )
        },
        "shortfall": mu.shortfall(a.get("blocks_tested_against_a_coding_gene") or 0, "gene"),
        "shortfall_upper_bound": mu.shortfall(b["blocks_touched"], "gene"),
    }


#: how the files a run reads under data/ are grouped in its manifest, first prefix that matches
READ_GROUPS = (
    ("data/knowledge/alphagenome/all_elements/", "the element archive the population is built from"),
    ("data/knowledge/measured_unknown/", "the public files this lane fetched, digested member by member"),
    ("data/results/", "the organiser's per-chromosome inputs, as committed"),
    ("data/reference/", "the reference genome and annotation (coding genes; the guide build check)"),
    ("data/", "the held assay tables (benchmark, lentiMPRA joint library, VISTA, growth-screen guides)"),
)


def result_manifest(sources: list[dict[str, Any]], reads: list[str]) -> dict[str, Any]:
    groups: dict[str, list[str]] = {}
    for path in sorted(set(reads)):
        label = next(lab for pre, lab in READ_GROUPS if path.startswith(pre))
        groups.setdefault(label, []).append(path)
    inputs = [mf.files_entry(label, paths, partition="read by the run") for label, paths in groups.items()]
    return {
        "sources": [
            {"accession": f"{s['key']}: {s['label']}", "version": (s["query"] or s["licence"])[:400]}
            for s in sources
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38; hg19 tables lifted with UCSC hg19ToHg38.over.chain at liftOver's default "
        "minMatch 0.95; a file with no stated build is decided by its guide sequences against GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "reciprocal_overlap": mu.RECIPROCAL_OVERLAP,
            "min_match": mu.MIN_MATCH,
            "floor_blocks": mu.FLOOR_BLOCKS,
            "sizing": {k: list(v) if v else None for k, v in mu.SIZING.items()},
            "min_log2": mu.MIN_LOG2,
            "assembly_check": {"sample": ASSEMBLY_SAMPLE, "agree": ASSEMBLY_AGREE},
        },
        "exclusions": [
            "every readout column of every source: effect sizes, activities, p-values, calls and PIPs",
            "the fine-mapped eQTL track itself (gtexDapg), which is the readout",
            "mouse (mm10) files, and unstated-build files whose guides match neither GRCh38 nor lifted hg19",
        ],
        "partitions": "n/a: an inventory, nothing is trained or held out",
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def table_row(s: dict[str, Any]) -> dict[str, Any]:
    c = s.get("census") or {}
    return {
        "source": s["key"],
        "endpoint": s["endpoint"],
        "held_before_this_lane": s["held_before_this_lane"],
        "blocks_touched_of_882": c.get("blocks_touched"),
        "blocks_touched_of_531_carrying": c.get("blocks_touched_carrying"),
        "blocks_element_matched_at_0_5": c.get("blocks_element_matched"),
        "blocks_tested_against_a_coding_gene": c.get("blocks_tested_against_a_coding_gene"),
        "blocks_tested_against_an_attributed_gene": c.get("blocks_tested_against_an_attributed_gene"),
        "megabytes": round(s["bytes"] / 1e6, 2),
        "licence": s["licence"],
        "not_computed": s["not_computed"] if s.get("census") is None else None,
    }


#: what a person saw of any readout while this lane was built, kept so the "0 values read" is checkable
INCIDENTAL_DISPLAYS = [
    "2026-10-10, a header peek printed two data rows of ENCFF528RAY (PPIF, hg19) and three of ENCFF149IDL "
    "(FADS2, GRCh38), effect-size columns included; a later header check printed one row of ENCFF197SZI "
    "(a MYC guide's count). All six intervals lie outside every one of the 882 blocks (checked on the "
    "lifted coordinates), and none enters any count",
]


def verdict(sources: list[dict[str, Any]], union: dict[str, Any]) -> dict[str, Any]:
    """The inventory's answer, written from its own counts with no word chosen after them."""
    a = union["gene_resolved_sources"]
    coding = a.get("blocks_tested_against_a_coding_gene") or 0
    upper = union["with_unresolved_libraries_as_upper_bound"]["blocks_touched"]
    gene_ok = coding >= mu.FLOOR_BLOCKS
    eqtl = next(s for s in sources if s["key"] == "gtex_v8_finemapped_eqtl")["census"]
    activity = {
        s["key"]: s["census"]["blocks_touched"]
        for s in sources
        if s["endpoint"] == "activity" and s.get("census")
    }
    cell = {
        s["key"]: s["census"]["blocks_touched"]
        for s in sources
        if s["endpoint"] == "cell" and s.get("census")
    }
    return {
        "element_perturbation_suffices": gene_ok,
        "element_perturbation": (
            "Every public element-perturbation screen with a gene readout, pooled, tested an element "
            f"against a protein-coding gene in {coding} of the 882 real-unknown blocks ({upper} counting "
            "the libraries whose gene is not resolved per guide, as an upper bound). The smallest "
            f"calibrated design needs {mu.FLOOR_BLOCKS} compared blocks, so the shortfall is "
            f"{max(0, mu.FLOOR_BLOCKS - coding)} blocks to that floor and {max(0, 531 - coding)} to the "
            "design for a half; it is " + ("met." if gene_ok else "not met by anything published.")
        ),
        "reporter_and_growth": (
            "Reporter assays reach more blocks ("
            + ", ".join(f"{k} {v}" for k, v in activity.items())
            + ") and the growth screens "
            + ", ".join(f"{k} {v}" for k, v in cell.items())
            + ", but neither names a gene, so neither can score the clause's 'attributed to a gene'."
        ),
        "natural_perturbation": (
            "The one source whose tested universe reaches the sizing is GTEx v8 fine-mapped cis-eQTL: "
            f"{eqtl['blocks_touched']} of 882 blocks ({eqtl['blocks_touched_carrying']} of the 531 "
            "carrying a scored element) hold at least one dbSNP common single-nucleotide variant, the "
            "proxy used here for a variant GTEx tested. Its endpoint names a gene and a tissue. Its "
            "sensitivity rises with a variant's frequency, which purifying selection lowers, so a fair "
            "test must match its controls on the variants they carry; the fine-mapped track itself has "
            "not been read."
        ),
        "decision": {
            "registered_test": False,
            "why": (
                "Among the sources inventoried here, only CRISPR perturbation names the gene an element "
                f"was tested against, and it reaches {coding} blocks against a coding gene ({upper} as an "
                f"upper bound), under the {mu.FLOOR_BLOCKS} the smallest calibrated design needs. The "
                "calibrated sizing (cab70d9) describes that endpoint and no other. The fine-mapped eQTL "
                f"route reaches {eqtl['blocks_touched']} blocks by the tested-variant proxy, but its "
                "window rate, its clustering and its sensitivity have never been calibrated in this "
                "project, and its control has yet to be chosen, so the sizing cannot say it is enough. "
                "No test is registered by this lane."
            ),
            "what_the_clause_needs": [
                f"{max(0, mu.FLOOR_BLOCKS - coding)} more blocks with an element tested against a coding "
                f"gene for the smallest calibrated design, {max(0, 531 - coding)} more for the design for a "
                "half; no published screen inventoried here supplies them",
                "or, on natural variation: a power calibration of a fine-mapped eQTL arm measured on "
                "sequence outside the 882 blocks, a control matched on length, common-variant count and "
                "distance to the named gene's start, and a registration, all before the gtexDapg track "
                "is read over the blocks",
            ],
            "prior_exposure": (
                "earlier lanes range-read the gtexDapg track over storage units of the UNKNOWN space "
                "(human_panel) and at 1,505 CRISPRi elements (finemap_dapg); whether any of those rows "
                "lie in the 882 blocks is not checked here, because checking would read them"
            ),
        },
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=(__doc__ or "").split("\n")[0])
    ap.add_argument("--no-save", action="store_true")
    ap.add_argument("--offline", action="store_true", help="use only files already cached")
    args = ap.parse_args(argv)
    mf.trace_begin()
    t0 = time.time()
    rows = mu.population()
    g = mu.gate(rows)
    say("population", g)
    if not g["ok"]:
        say("the population is not the clause-2 population; nothing is counted")
        return 2
    coding = coding_symbols()
    say("coding symbols", len(coding))
    chain_file = fetch(CHAIN_URL, "hg19ToHg38.over.chain.gz", None, args.offline)
    chain = mu.Chain.from_file(Path(chain_file["path"]))
    sources: list[dict[str, Any]] = [held_benchmark(*b) for b in BENCHMARK]
    sources += flowfish(chain, args.offline)
    sources.append(gasperini(chain, args.offline))
    sources.append(scrna_guides(chain, args.offline))
    sources.append(facs(args.offline))
    sources.append(cell_level(args.offline))
    sources += mpra_sources(chain, args.offline)
    sources.append(vista_source())
    sources.append(mpravardb(rows))
    sources.append(gtex_tested_universe(rows))
    sources += not_fetched()
    for s in sources:
        if s["_rows"] is not None:
            s["census"] = mu.census(rows, s["_rows"], coding if s["gene_resolved_per_interval"] else None)
            basis = s["census"].get("blocks_tested_against_a_coding_gene")
            s["shortfall"] = mu.shortfall(
                s["census"]["blocks_touched"] if basis is None else basis, s["endpoint"]
            )
        say(s["key"], table_row(s))
    union = gene_union(rows, sources, coding)
    say(
        "gene-level union",
        {k: v for k, v in union["gene_resolved_sources"].items() if not isinstance(v, list)},
    )
    say(
        "upper bound",
        {
            k: v
            for k, v in union["with_unresolved_libraries_as_upper_bound"].items()
            if not isinstance(v, list)
        },
    )
    for s in sources:
        s.pop("_rows", None)
    record = {
        "result": NAME,
        "date": time.strftime("%Y-%m-%d"),
        "question": "does any public measurement reach enough of the 882 real-unknown blocks to score "
        "milestone 1.3's clause 2 (a gene and a tissue per block) at the calibrated sizing?",
        "population": {
            "gate": g,
            "by_case": dict(Counter(r["case"] for r in rows)),
            "blocks_naming_a_coding_gene": sum(1 for r in rows if r["named_coding_genes"]),
            "blocks_naming_any_gene": sum(1 for r in rows if r["named_genes"]),
        },
        "sizing": {
            "source": "lane-s2 cab70d9 (registered 019a9ab), calibrated on CRISPRi element-gene readouts",
            "designs": {k: list(v) if v else None for k, v in mu.SIZING.items()},
            "floor_blocks": mu.FLOOR_BLOCKS,
        },
        "readout_values_read": 0,
        "incidental_displays": INCIDENTAL_DISPLAYS,
        "verdict": verdict(sources, union),
        "table": [table_row(s) for s in sources],
        "gene_level_union": union,
        "sources": sources,
        "lift": {"chain": chain_file},
        "alphagenome_requests": 0,
        "money": "none",
        "seconds": round(time.time() - t0, 1),
    }
    record["result_manifest"] = result_manifest(sources, mf.traced_reads())
    if args.no_save:
        json.dump(
            {k: record[k] for k in ("population", "verdict", "table")}, sys.stdout, indent=1, default=str
        )
        return 0
    path = save_result(NAME, record)
    say("wrote", path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
