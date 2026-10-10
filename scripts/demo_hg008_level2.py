# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""A demo tumour at data level 2: GIAB HG008 with its matched normal genome.

The demo tumour (`data/demo/cancer_tumour.vcf`) runs at data level 1. Level 2 is
the matched normal's own variant file, held by the profile and subtracted by
`rank_variants`; a somatic call set made against a normal elsewhere is level-1
data with a provenance claim. TCGA and CPTAC publish no normal genome openly, so
neither can reach level 2. docs/THERAPEUTICS.md, "A demo case above data level
1", records that decision and this case, both committed before this script was
written.

The case: GIAB HG008, a pancreatic ductal adenocarcinoma cell line (HG008-T,
Liss lab batch 0823p23) and the same patient's normal duodenal tissue
(HG008-N-D), consented for fully public release. The inputs are two files of
the NYGC Somatic Pipeline v7 analysis (2024-04-12, GRCh38-GIABv3):

  tumour  HG008-T--HG008-N.snv.indel.high_confidence.v7.annotated.vcf.gz
          the HG008-T column of NYGC's high-confidence somatic calls
  normal  HG008-N.haplotypecaller.gatk.annotated.unfiltered.vcf.gz
          every HaplotypeCaller call in HG008-N, common variants included

Both are checked against NIST's published md5 and cut to the same region set,
the union of every GENCODE v50 CDS feature padded by 2 bp on each side; a record
is inside when any base of its reference allele is. The 1.5 GB normal file is
streamed and never stored: it is hashed and filtered as it arrives.

Two commands:

  build  writes data/demo/hg008_pdac/HG008-T.vcf, HG008-N-D.vcf and
         provenance.json (sources, md5s, the region rule and its size, what
         each rewrite dropped).
  run    runs the therapeutic pipeline on the pair (level 2) and on the tumour
         file alone (level 1), writes the level-2 analysis to
         data/demo/hg008_pdac/therapeutics/ and the comparison to
         data/demo/hg008_pdac/level_comparison.json. New VEP annotations go to
         data/demo/hg008_pdac/vep_cache.jsonl, so a later run reads them from
         the demo itself rather than from the network.

Run: `uv run python scripts/demo_hg008_level2.py build`
     `uv run python scripts/demo_hg008_level2.py run [--offline]`
"""

from __future__ import annotations

import bisect
import gzip
import hashlib
import io
import json
import sys
import urllib.request
from datetime import date
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

DEMO = ROOT / "data/demo/hg008_pdac"
TUMOUR = DEMO / "HG008-T.vcf"
NORMAL = DEMO / "HG008-N-D.vcf"
PROVENANCE = DEMO / "provenance.json"
COMPARISON = DEMO / "level_comparison.json"
OUTPUTS = DEMO / "therapeutics"

BASE = (
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/analysis/"
    "NYGC-somatic-pipeline_20240412/GRCh38-GIABv3/"
)
SOURCES = {
    "tumour": {
        "file": "HG008-T--HG008-N.snv.indel.high_confidence.v7.annotated.vcf.gz",
        "md5": "72e85c9ddd6c222137e4136f6a76d21f",
        "sample": "HG008-T",
    },
    "normal": {
        "file": "HG008-N.haplotypecaller.gatk.annotated.unfiltered.vcf.gz",
        "md5": "891836dd9da087bceea43ee3e436f1dc",
        "sample": "HG008-N",
    },
}
GENCODE = ROOT / "data/reference"
CHROMS = [f"chr{c}" for c in [*range(1, 23), "X", "Y", "M"]]
PAD = 2
LICENCE = (
    "NIST data use policy: not subject to copyright in the United States (17 U.S.C. 105); "
    "royalty-free worldwide grant elsewhere. Modified works must state the change and acknowledge "
    "the National Institute of Standards and Technology as the source of the data."
)
CHANGE_NOTICE = (
    "Changed by GenomeOS on {day}: cut to GENCODE v50 CDS +-2 bp, reduced to one sample, INFO dropped "
    "(scripts/demo_hg008_level2.py). Source: GIAB HG008, NYGC Somatic Pipeline v7, NIST."
)


# ---- the region set -----------------------------------------------------------------


def region_set() -> tuple[dict[str, tuple[list[int], list[int]]], dict[str, Any]]:
    """Every GENCODE v50 CDS feature, padded by PAD bp and merged, per chromosome (1-based, closed)."""
    merged: dict[str, tuple[list[int], list[int]]] = {}
    files: dict[str, str] = {}
    total = 0
    for chrom in CHROMS:
        path = GENCODE / f"gencode_v50_{chrom}.gff3.gz"
        if not path.exists():
            raise SystemExit(f"missing {path.relative_to(ROOT)}: the region set needs every chromosome")
        files[str(path.relative_to(ROOT))] = "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()
        spans: list[tuple[int, int]] = []
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split("\t", 5)
                if len(f) > 4 and f[2] == "CDS":
                    spans.append((max(1, int(f[3]) - PAD), int(f[4]) + PAD))
        spans.sort()
        starts: list[int] = []
        ends: list[int] = []
        for s, e in spans:
            if starts and s <= ends[-1] + 1:
                ends[-1] = max(ends[-1], e)
            else:
                starts.append(s)
                ends.append(e)
        merged[chrom] = (starts, ends)
        total += sum(e - s + 1 for s, e in zip(starts, ends, strict=True))
    meta = {
        "rule": (
            f"union of every GENCODE v50 CDS feature, padded by {PAD} bp on each side and merged; a VCF "
            "record is inside when any base of its reference allele is"
        ),
        "chromosomes": CHROMS,
        "intervals": sum(len(s) for s, _ in merged.values()),
        "bases": total,
        "gencode_files": files,
    }
    return merged, meta


def inside(regions: dict[str, tuple[list[int], list[int]]], chrom: str, pos: int, ref: str) -> bool:
    if chrom not in regions:
        return False
    starts, ends = regions[chrom]
    last = pos + max(len(ref), 1) - 1
    i = bisect.bisect_right(starts, last) - 1
    return i >= 0 and ends[i] >= pos


# ---- streaming with a checksum ---------------------------------------------------------


class _Hashing(io.RawIOBase):
    """A readable stream that hashes every byte it passes on."""

    def __init__(self, raw: Any) -> None:
        self.raw = raw
        self.md5 = hashlib.md5()  # noqa: S324 - NIST publishes md5; this checks it, it secures nothing
        self.bytes = 0

    def readable(self) -> bool:
        return True

    def readinto(self, buffer: Any) -> int:
        data = self.raw.read(len(buffer))
        n = len(data)
        buffer[:n] = data
        self.md5.update(data)
        self.bytes += n
        return n

    def drain(self) -> None:
        while self.readinto(bytearray(1 << 20)):
            pass


def _stream(role: str):
    src = SOURCES[role]
    resp = urllib.request.urlopen(BASE + src["file"], timeout=600)  # noqa: S310 - fixed https source
    hashing = _Hashing(resp)
    text = io.TextIOWrapper(gzip.GzipFile(fileobj=io.BufferedReader(hashing, 1 << 20)), encoding="utf-8")
    return hashing, text


def _verify(role: str, hashing: _Hashing) -> dict[str, Any]:
    hashing.drain()
    got = hashing.md5.hexdigest()
    want = SOURCES[role]["md5"]
    if got != want:
        raise SystemExit(f"{role}: md5 {got} does not match NIST's {want}; nothing written")
    return {"url": BASE + SOURCES[role]["file"], "md5": got, "bytes": hashing.bytes}


def _header(sample: str, fmt_lines: list[str], role: str, day: str) -> list[str]:
    return [
        "##fileformat=VCFv4.2",
        f"##source=GIAB HG008 {role}, NYGC Somatic Pipeline v7 (2024-04-12), GRCh38-GIABv3",
        "##reference=GRCh38",
        f"##GenomeOS_change={CHANGE_NOTICE.format(day=day)}",
        f"##GenomeOS_licence={LICENCE}",
        *fmt_lines,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample,
    ]


# ---- build -----------------------------------------------------------------------------


def _tumour(regions, day: str) -> dict[str, Any]:
    hashing, text = _stream("tumour")
    col = None
    kept: list[str] = []
    seen = outside = 0
    af_gap = 0.0
    for line in text:
        if line.startswith("##"):
            continue
        f = line.rstrip("\n").split("\t")
        if line.startswith("#CHROM"):
            col = f.index(SOURCES["tumour"]["sample"])
            continue
        seen += 1
        if not inside(regions, f[0], int(f[1]), f[3]):
            outside += 1
            continue
        sample = dict(zip(f[8].split(":"), f[col].split(":"), strict=False))
        ad = [int(x) for x in sample["AD"].split(",")]
        af = ad[1] / sum(ad) if sum(ad) else 0.0
        if sample.get("AF") not in (None, "."):
            af_gap = max(af_gap, abs(af - float(sample["AF"].split(",")[0])))
        values = f"{sample['AD']}:{sample['DP']}:{af:.4f}"
        kept.append("\t".join([f[0], f[1], ".", f[3], f[4], ".", f[6], ".", "AD:DP:AF", values]))
    source = _verify("tumour", hashing)
    fmt = [
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths in HG008-T, as NYGC reports">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth in HG008-T, as NYGC reports it">',
        '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Alternate allele fraction, alt/(ref+alt) of AD">',
    ]
    role = "tumour cell line HG008-T (batch 0823p23), somatic calls against HG008-N-D"
    TUMOUR.write_text("\n".join([*_header("HG008-T", fmt, role, day), *kept]) + "\n")
    return {
        **source,
        "records_in_source": seen,
        "records_outside_region_set": outside,
        "records_written": len(kept),
        "genotype": "the source records carry no GT; none is written",
        "allele_fraction": "alt / (ref + alt) of the tumour's AD, rounded to 4 places",
        "largest_difference_from_source_af": round(af_gap, 6),
    }


def _normal(regions, day: str) -> dict[str, Any]:
    hashing, text = _stream("normal")
    kept: list[str] = []
    seen = outside = no_alt = multi = 0
    filters: dict[str, int] = {}
    for line in text:
        if line.startswith("#"):
            continue
        seen += 1
        head = line.split("\t", 7)
        chrom, pos, ref, alt, flt = head[0], int(head[1]), head[3], head[4], head[6]
        if not inside(regions, chrom, pos, ref):
            outside += 1
            continue
        fmt, sample = line.rstrip("\n").rsplit("\t", 2)[1:]
        keys = fmt.split(":")
        gt = sample.split(":")[keys.index("GT")] if "GT" in keys else "."
        alleles = [a for a in gt.replace("|", "/").split("/") if a not in (".", "")]
        if not any(a != "0" for a in alleles):
            no_alt += 1
            continue
        if "," in alt:
            multi += 1
        filters[flt] = filters.get(flt, 0) + 1
        kept.append("\t".join([chrom, str(pos), ".", ref, alt, ".", flt, ".", "GT", gt]))
    source = _verify("normal", hashing)
    fmt_lines = [
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype in HG008-N, as HaplotypeCaller called">'
    ]
    role = "matched normal HG008-N-D (duodenal tissue), all germline calls"
    NORMAL.write_text("\n".join([*_header("HG008-N-D", fmt_lines, role, day), *kept]) + "\n")
    return {
        **source,
        "records_in_source": seen,
        "records_outside_region_set": outside,
        "records_without_an_alternate_allele_dropped": no_alt,
        "records_written": len(kept),
        "multiallelic_records_written": multi,
        "filter_values_written": dict(sorted(filters.items(), key=lambda kv: -kv[1])),
        "note": (
            "FILTER is kept as HaplotypeCaller wrote it; the pipeline's subtraction reads every record "
            "whatever its FILTER, and matches chromosome, position, reference and the whole alternate list"
        ),
    }


def build() -> int:
    DEMO.mkdir(parents=True, exist_ok=True)
    day = date.today().isoformat()
    regions, meta = region_set()
    print(f"region set: {meta['intervals']:,} intervals, {meta['bases']:,} bases", flush=True)
    tumour = _tumour(regions, day)
    print(f"tumour: {tumour['records_written']} of {tumour['records_in_source']} records inside", flush=True)
    normal = _normal(regions, day)
    print(f"normal: {normal['records_written']} of {normal['records_in_source']} records inside", flush=True)
    provenance = {
        "case": (
            "GIAB HG008: pancreatic ductal adenocarcinoma cell line HG008-T (Liss lab batch 0823p23, "
            "passage 23) "
            "and matched normal duodenal tissue HG008-N-D, one female patient"
        ),
        "consent": (
            "NIST: 'explicitly consented for fully public dissemination of genomic data and cell lines'"
        ),
        "analysis": (
            "NYGC Somatic Pipeline v7, 2024-04-12, GRCh38-GIABv3 (README_NYGC.md in the source directory)"
        ),
        "licence": LICENCE,
        "built": day,
        "script": "scripts/demo_hg008_level2.py build",
        "decided_in": (
            "docs/THERAPEUTICS.md, 'A demo case above data level 1, decided 2026-10-10 before it was built'"
        ),
        "region_set": meta,
        "tumour": {"written_to": str(TUMOUR.relative_to(ROOT)), **tumour},
        "normal": {"written_to": str(NORMAL.relative_to(ROOT)), **normal},
        "not_supplied": {
            "tumour_rna": "raw FASTQ only (UMD 2025-09-25), from NIST passages 21 and 103, not batch 0823p23",
            "copy_number": "NYGC BIC-seq2 segments exist; not supplied",
            "hla": "Kourami alleles on HG008-N exist; not supplied",
        },
    }
    PROVENANCE.write_text(json.dumps(provenance, indent=1) + "\n")
    return 0


# ---- run -------------------------------------------------------------------------------


def _counts(tumour: str, normal: str | None) -> dict[str, int]:
    from genomeos.cancer.compare import somatic
    from genomeos.genome.variants import iter_vcf

    total = sum(1 for _ in iter_vcf(tumour, None, pass_only=False))
    kept = len(somatic(normal, tumour)) if normal else total
    return {"tumour_records": total, "removed_by_the_normal": total - kept, "records_ranked": kept}


def _summary(a: dict[str, Any], tumour: str, normal: str | None) -> dict[str, Any]:
    level = a["data_level"]
    return {
        "level_reached": level["level_reached"],
        "highest_input_present": level["highest_input_present"],
        "inputs_present": level["inputs_present"],
        **_counts(tumour, normal),
        "variants_total": a["variants_total"],
        "somatic_candidates": a["somatic_candidates"],
        "likely_germline": a["variants_total"] - a["somatic_candidates"],
        "missing_data": [m["input"] for m in a["missing_data"]],
        "candidates": [
            {
                "gene": c.gene,
                "target_class": c.target_class,
                "evidence_tier": c.evidence_tier,
                "somatic_status": [o.somatic_status for o in c.origins],
                "expression_known": c.tumour.expression.known,
            }
            for c in a["candidates"]
        ],
    }


def run(net: bool) -> int:
    from genomeos.therapeutics import analyse_vcf, write_outputs
    from genomeos.therapeutics.report import global_limitations

    if not (TUMOUR.exists() and NORMAL.exists()):
        raise SystemExit("run `build` first")
    tumour, normal = str(TUMOUR.relative_to(ROOT)), str(NORMAL.relative_to(ROOT))
    common = {"net": net, "indirect": True, "vep_append_to": DEMO, "log": sys.stderr}
    level2 = analyse_vcf(tumour, normal_vcf=normal, **common)
    level1 = analyse_vcf(tumour, **common)
    write_outputs(level2, OUTPUTS)
    out = {
        "case": "GIAB HG008 (data/demo/hg008_pdac/provenance.json)",
        "mode": "network" if net else "offline",
        "level_2": _summary(level2, tumour, normal),
        "level_1_same_tumour_file": _summary(level1, tumour, None),
        "refused_at_level_2": global_limitations(level2),
    }
    COMPARISON.write_text(json.dumps(out, indent=1, default=str) + "\n")
    for k in ("level_2", "level_1_same_tumour_file"):
        s = out[k]
        print(
            f"{k}: level {s['level_reached']}, {s['tumour_records']} records, "
            f"{s['removed_by_the_normal']} removed by the normal, {s['likely_germline']} set aside as "
            f"likely germline, {len(s['candidates'])} candidates",
            flush=True,
        )
    return 0


def main(argv: list[str]) -> int:
    if not argv or argv[0] not in {"build", "run"}:
        print(__doc__)
        return 2
    if argv[0] == "build":
        return build()
    return run(net="--offline" not in argv)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
