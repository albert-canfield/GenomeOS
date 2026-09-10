"""VCF reading and variant application (task 1.3).

Streams a VCF (plain or gzip), keeps the first sample's genotype, and applies
the chosen haplotype's alleles to a reference sequence to produce an
individual's chromosome. Phased genotypes (0|1) keep their haplotype order;
unphased heterozygous calls (0/1) put the ALT on haplotype 1 by default.
"""

from __future__ import annotations

import gzip
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path

from .sequence import Sequence


@dataclass(frozen=True, slots=True)
class Variant:
    chrom: str
    pos: int  # 0-based
    ref: str
    alts: tuple[str, ...]
    id: str = "."
    gt: tuple[int, ...] = ()  # allele indexes per haplotype; () if missing
    phased: bool = False
    info: str = ""

    @property
    def is_snv(self) -> bool:
        return len(self.ref) == 1 and all(len(a) == 1 for a in self.alts)

    def allele(self, haplotype: int) -> str | None:
        """Allele carried on haplotype 0 or 1; None if reference or missing."""
        if len(self.gt) <= haplotype:
            return None
        idx = self.gt[haplotype]
        if idx <= 0:
            return None
        return self.alts[idx - 1]


def _open(path: str | Path):
    path = Path(path)
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)


def iter_vcf(path: str | Path, chroms: set[str] | None = None, pass_only: bool = True) -> Iterator[Variant]:
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if chroms is not None and f[0] not in chroms:
                continue
            if pass_only and f[6] not in ("PASS", "."):
                continue
            gt: tuple[int, ...] = ()
            phased = False
            if len(f) >= 10:
                keys = f[8].split(":")
                vals = f[9].split(":")
                if "GT" in keys:
                    raw = vals[keys.index("GT")]
                    phased = "|" in raw
                    parts = raw.replace("|", "/").split("/")
                    if all(p not in (".", "") for p in parts):
                        gt = tuple(int(p) for p in parts)
            yield Variant(f[0], int(f[1]) - 1, f[3], tuple(f[4].split(",")), f[2], gt, phased, f[7])


def apply_variants(
    reference: Sequence, variants: Iterable[Variant], haplotype: int = 0
) -> tuple[Sequence, dict]:
    """Return the haplotype sequence and statistics. Overlapping variants after
    the first are skipped; reference mismatches are skipped and counted."""
    ref = str(reference)
    chunks: list[str] = []
    cursor = 0
    stats = {
        "applied": 0,
        "snv": 0,
        "indel": 0,
        "skipped_overlap": 0,
        "skipped_ref_mismatch": 0,
        "length_delta": 0,
    }
    for v in sorted(variants, key=lambda x: x.pos):
        alt = v.allele(haplotype)
        if alt is None:
            continue
        if v.pos < cursor:
            stats["skipped_overlap"] += 1
            continue
        if ref[v.pos : v.pos + len(v.ref)].upper() != v.ref.upper():
            stats["skipped_ref_mismatch"] += 1
            continue
        chunks.append(ref[cursor : v.pos])
        chunks.append(alt)
        cursor = v.pos + len(v.ref)
        stats["applied"] += 1
        stats["snv" if len(v.ref) == 1 and len(alt) == 1 else "indel"] += 1
        stats["length_delta"] += len(alt) - len(v.ref)
    chunks.append(ref[cursor:])
    return Sequence("".join(chunks)), stats


def write_haplotypes(
    reference: Sequence, chrom: str, variants: list[Variant], out_fasta: Path, sample: str
) -> dict:
    """Write hap1/hap2 records to a FASTA and return per-haplotype stats."""
    result = {}
    with open(out_fasta, "w") as fh:
        for h in (0, 1):
            seq, stats = apply_variants(reference, variants, h)
            result[f"hap{h + 1}"] = stats
            header = f">{chrom}_{sample}_hap{h + 1} {chrom} haplotype {h + 1} of {sample}"
            fh.write(f"{header}, {stats['applied']} variants\n")
            s = str(seq)
            for i in range(0, len(s), 60):
                fh.write(s[i : i + 60] + "\n")
    return result
