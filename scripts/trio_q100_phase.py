# SPDX-License-Identifier: AGPL-3.0-or-later
"""Place HG002's trio de novo candidates on the phased, parent-labelled Q100 assembly (docs/DATA.md,
"Phasing the trio candidates on the Q100 assembly", pre-registered 2026-09-27 before this ran).

Inputs are local and git-ignored: the three individuals under data/individuals (HG002, HG003, HG004
with their trusted regions) and NIST/GIAB's v5.0q release files under data/cache/q100, fetched from
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/
(dipcall_output/GRCh38_HG2-T2TQ100-V1.1_dipcall-z2k.dip.vcf.gz and .dip.bed, and
HG002_GRCh38_v5.0q_smvar.benchmark.bed). Data: NIST, not subject to US copyright (17 USC 105); the
source is acknowledged here and in docs/DATA.md.

Writes data/results/trio_q100_phase.json: counts per class, never a position, allele or genotype.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.genome import individuals as ind  # noqa: E402
from genomeos.results import save_result  # noqa: E402

CACHE = Path("data/cache/q100")
VCF = CACHE / "GRCh38_HG2-T2TQ100-V1.1_dipcall-z2k.dip.vcf.gz"
DIP = CACHE / "GRCh38_HG2-T2TQ100-V1.1_dipcall-z2k.dip.bed"
BENCH = CACHE / "HG002_GRCh38_v5.0q_smvar.benchmark.bed"
# NIST's HG002 de novo and mosaic regions excluded from v5.0q, named in defrabb_files/resources.yml
MOSAIC = CACHE / "mergedQ100denovomosaicexclusions_240506_GRCh38_withoverlaprepeats.bed.gz"
OUT = Path("data/results/trio_q100_phase.json")  # written through save_result (item 12 S6)
REGISTERED_CANDIDATES = 1430  # the trio's normalised count, docs/DATA.md 2026-09-12
TRIO = ("HG002", "HG003", "HG004")
AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
V5 = (
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/"
)
SOURCES = [
    {
        "accession": "GIAB AshkenazimTrio HG002 v5.0q: dipcall_output/"
        "GRCh38_HG2-T2TQ100-V1.1_dipcall-z2k.dip.vcf.gz and .dip.bed, HG002_GRCh38_v5.0q_smvar.benchmark.bed",
        "version": "v5.0q (T2T HG002 Q100 v1.1), md5-checked against checksum.md5",
        "url": V5,
        "licence": "NIST, not subject to US copyright (17 USC 105), AS IS",
    },
    {
        "accession": "NIST mergedQ100denovomosaicexclusions_240506_GRCh38_withoverlaprepeats.bed.gz",
        "version": "240506, named in the v5.0q defrabb_files/resources.yml",
        "licence": "NIST, not subject to US copyright (17 USC 105), AS IS",
    },
    {
        "accession": "GIAB v4.2.1 benchmark calls and regions for HG002, HG003 and HG004 (GRCh38)",
        "version": "NISTv4.2.1",
        "licence": "NIST, not subject to US copyright (17 USC 105), AS IS",
    },
    {"accession": "UCSC hg38 chromosome FASTA (GRCh38), data/reference", "version": "hg38"},
]


def trio_inputs(chroms: list[str] = AUTOSOMES) -> list[dict]:
    """The trio's own calls and trusted regions, and the reference bases read, one entry per person."""
    from genomeos.genome import reference_fasta

    out = []
    for name in TRIO:
        files = [ind.vcf_path(name, c) for c in chroms]
        files += [ind.ROOT / name / f"regions_{c}.bed" for c in chroms]
        out.append(
            mf.files_entry(
                f"{name}: GIAB v4.2.1 calls and trusted regions, {chroms[0]}-{chroms[-1]}",
                [f for f in files if f is not None and Path(f).exists()],
            )
        )
    out.append(
        mf.files_entry(
            f"data/reference GRCh38 {chroms[0]}-{chroms[-1]}",
            [reference_fasta(c, Path("data/reference")) for c in chroms],
        )
    )
    return out


def manifest(parameters: dict, extra_inputs: list[dict] = (), extra_sources: list[dict] = ()) -> dict:
    return {
        "sources": SOURCES + list(extra_sources),
        "inputs": [mf.input_entry(p) for p in (VCF, DIP, BENCH, MOSAIC)] + trio_inputs() + list(extra_inputs),
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": parameters,
        "exclusions": [
            "calls outside the trusted regions of all three people or outside the dipcall .dip.bed",
            "NIST's de novo and mosaic exclusion regions, counted apart rather than dropped",
        ],
        "partitions": "n/a: two measurements of one person compared; no evaluation split",
    }


def main() -> int:
    t0 = time.time()
    r = ind.phase_trio_by_assembly(
        "HG002",
        "HG003",
        "HG004",
        VCF,
        DIP,
        BENCH,
        hap1_parent="father",
        other_beds={"nist_de_novo_mosaic_exclusions": MOSAIC},
    )
    if r["candidates"]["total"] != REGISTERED_CANDIDATES:
        print(f"refused: {r['candidates']['total']} candidates, registered {REGISTERED_CANDIDATES}")
        return 1
    r["source"] = {
        "assembly": "T2T HG002 Q100 v1.1 (trio-binned, maternal and paternal haplotypes), dipcall against "
        "GRCh38 as released by NIST/GIAB in v5.0q (dipcall_output/, md5-checked)",
        "haplotype_order": "defrabb passes paternal.fa as dipcall's first haplotype; checked on the file: "
        "every PASS chrX non-PAR call is .|1 and every chrY call 1|.",
        "licence": "NIST data, not subject to copyright in the US (17 USC 105), provided AS IS; "
        "acknowledgement: National Institute of Standards and Technology, Genome in a Bottle",
        "trio": "GIAB v4.2.1 benchmark calls and regions for HG002, HG003, HG004, alleles normalised",
    }
    r["registration"] = (
        "docs/DATA.md, section 'Phasing the trio candidates on the Q100 assembly (2026-09-27)'"
    )
    r["evidence"] = "measured: two independent measurements of one person compared (read-based benchmark "
    r["evidence"] += "calls against a phased assembly); counts only"
    r["seconds"] = round(time.time() - t0, 1)
    params = {
        "registered_candidates": REGISTERED_CANDIDATES,
        "hap1_parent": "father",
        "chromosomes": AUTOSOMES,
    }
    save_result(OUT.stem, r, manifest=manifest(params))
    print(
        json.dumps(
            {
                k: r[k]
                for k in (
                    "candidates",
                    "control",
                    "control_found_on_one_haplotype",
                    "control_parent_agreement",
                    "candidates_paternal_share",
                )
            },
            indent=1,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
