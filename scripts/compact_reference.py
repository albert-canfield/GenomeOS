#!/usr/bin/env python3
"""Keep each chromosome once, as blocked gzip, instead of twice.

The reference cache holds every chromosome as a flat `chrN.fa` for random
access and as a plain `chrN.fa.gz` for streaming, because plain gzip cannot be
seeked. Blocked gzip can (`genomeos/genome/bgzf.py`), so one file does both
jobs and the cache shrinks by about four fifths with nothing lost.

Nothing is deleted on trust. For each chromosome this writes the blocked file
beside the old ones, checks that it decompresses to the SHA-256 of the flat
`.fa` **and** to the SHA-256 of the plain `.fa.gz` -- so both files that are
about to go are proved redundant, not assumed to be -- records those hashes in
`data/results/reference_bgzf.json`, and only then removes them. A chromosome a
process currently holds open is skipped.

    scripts/compact_reference.py --list
    scripts/compact_reference.py --chrom chr21 --dry-run
    scripts/compact_reference.py --all --skip chr2
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from genomeos.genome import bgzf  # noqa: E402
from genomeos.genome.index import write_fai  # noqa: E402

REFERENCE = Path("data/reference")
MANIFEST = Path("data/results/reference_bgzf.json")
LEVEL = 9  # level 6 costs 2.3% more disk for half the time; the file is written once


def chromosomes() -> list[str]:
    """Every chromosome with a flat `.fa` that could be folded into one file."""
    return sorted(p.name[:-3] for p in REFERENCE.glob("*.fa"))


def held_open() -> set[Path]:
    """Reference files some running process has open, which must not be touched.

    A chromosome being read right now (a scoring run, the web server) keeps its
    files: a reader that opened the flat `.fa` before the swap would lose it.
    """
    files = [p.resolve() for p in REFERENCE.iterdir() if p.suffix in (".fa", ".gz")]
    if not files:
        return set()
    try:
        out = subprocess.run(
            ["lsof", "-F", "n", "--", *[str(p) for p in files]],
            capture_output=True,
            text=True,
            timeout=300,
        ).stdout
    except (OSError, subprocess.SubprocessError) as ex:
        raise SystemExit(f"cannot ask lsof what is open ({ex}); refusing to delete anything") from ex
    return {Path(line[1:]) for line in out.splitlines() if line.startswith("n/")}


def sha256_of_gzip(path: Path) -> str:
    """SHA-256 of what a gzip file decompresses to, streamed."""
    h = hashlib.sha256()
    with gzip.open(path, "rb") as fh:
        while True:
            b = fh.read(1 << 22)
            if not b:
                return h.hexdigest()
            h.update(b)


def load_manifest() -> dict:
    if MANIFEST.exists():
        return json.loads(MANIFEST.read_text())
    return {
        "what": "Reference chromosomes folded from a flat .fa plus a plain .fa.gz into one blocked-gzip "
        ".fa.gz with a .gzi. The hashes are of the decompressed bytes, so tests/test_bgzf.py can "
        "re-prove that nothing was lost long after the originals were deleted.",
        "chromosomes": {},
    }


def save_manifest(m: dict) -> None:
    """Merge into whatever is on disk rather than overwrite it.

    Two runs can be in flight at once (the chromosomes take half an hour, the
    variant files a minute), and each holds a copy of the manifest it loaded at
    the start. Writing that copy back would drop the other run's rows -- which
    is exactly how the first pass lost the variant-file section.
    """
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    on_disk = json.loads(MANIFEST.read_text()) if MANIFEST.exists() else {}
    for key, value in m.items():
        if isinstance(value, dict) and isinstance(on_disk.get(key), dict):
            on_disk[key].update(value)
        else:
            on_disk[key] = value
    MANIFEST.write_text(json.dumps(on_disk, indent=1, sort_keys=True) + "\n")


def compact(chrom: str, dry_run: bool = False, keep: bool = False) -> dict | None:
    """Fold one chromosome. Returns its manifest row, or None when it was skipped."""
    flat = REFERENCE / f"{chrom}.fa"
    gz = REFERENCE / f"{chrom}.fa.gz"
    if not flat.exists() and gz.exists() and bgzf.is_bgzf(gz):
        print(f"{chrom}: already blocked")
        return None
    if not flat.exists():
        print(f"{chrom}: no flat .fa")
        return None
    open_now = held_open()
    if flat.resolve() in open_now or gz.resolve() in open_now:
        print(f"{chrom}: SKIPPED, a running process has it open")
        return None

    before = flat.stat().st_size + (gz.stat().st_size if gz.exists() else 0)
    print(f"{chrom}: {before / 1e6:.0f} MB in two files", flush=True)
    if dry_run:
        return None

    t0 = time.time()
    new = REFERENCE / f"{chrom}.fa.bgzf"
    bgzf.compress_file(flat, new, level=LEVEL)
    made = time.time() - t0

    flat_hash = bgzf.sha256_of(flat)
    new_hash = bgzf.sha256_of_bgzf(new)
    if new_hash != flat_hash:
        new.unlink(missing_ok=True)
        bgzf.gzi_path(new).unlink(missing_ok=True)
        raise SystemExit(f"{chrom}: blocked file does not match the .fa ({new_hash} != {flat_hash})")
    gz_hash = sha256_of_gzip(gz) if gz.exists() else flat_hash
    if gz_hash != flat_hash:
        new.unlink(missing_ok=True)
        bgzf.gzi_path(new).unlink(missing_ok=True)
        raise SystemExit(
            f"{chrom}: the flat .fa and the plain .fa.gz were never the same bytes "
            f"({flat_hash} != {gz_hash}); nothing deleted, this needs a person"
        )

    row = {
        "sha256": flat_hash,
        "uncompressed_bytes": flat.stat().st_size,
        "was_bytes": before,
        "now_bytes": new.stat().st_size + bgzf.gzi_path(new).stat().st_size,
        "level": LEVEL,
        "seconds": round(made, 1),
        "verified": "decompressed bytes of the blocked file match both the .fa and the plain .fa.gz",
    }

    # Replace in place: an open reader keeps its file, a new one finds the blocked file.
    fai = REFERENCE / f"{chrom}.fa.fai"
    if fai.exists():
        (REFERENCE / f"{chrom}.fa.gz.fai").write_bytes(fai.read_bytes())
    bgzf.gzi_path(new).rename(REFERENCE / f"{chrom}.fa.gz.gzi")
    new.rename(gz)
    if not (REFERENCE / f"{chrom}.fa.gz.fai").exists():
        write_fai(gz)
    if not keep:
        flat.unlink()
    saved = before - row["now_bytes"]
    print(
        f"{chrom}: {row['now_bytes'] / 1e6:.0f} MB in one file, {saved / 1e6:.0f} MB reclaimed "
        f"({made:.0f} s, sha256 {flat_hash[:12]})",
        flush=True,
    )
    return row


def compact_vcf(path: Path, keep: bool = False) -> dict | None:
    """Compress one split variant file, proving first that the same variants come back out.

    The proof is stronger than a byte hash alone: the file is read twice through
    `iter_vcf`, once from the plain file and once from the compressed one, and
    every variant must match in the same order, so a difference in what the
    readers *parse* is caught, not only a difference in bytes.
    """
    from genomeos.genome import iter_vcf

    gz = path.with_name(path.name + ".gz")
    if not path.exists():
        return None
    if path.resolve() in held_open():
        print(f"{path.name}: SKIPPED, a running process has it open")
        return None

    before = path.stat().st_size
    tmp = path.with_name(path.name + ".compacting.gz")  # ends in .gz so every reader opens it as one
    bgzf.compress_file(path, tmp, level=LEVEL)
    plain_hash = bgzf.sha256_of(path)
    if bgzf.sha256_of_bgzf(tmp) != plain_hash:
        tmp.unlink(missing_ok=True)
        bgzf.gzi_path(tmp).unlink(missing_ok=True)
        raise SystemExit(f"{path.name}: compressed bytes differ from the original")

    n = 0
    for a, b in zip(iter_vcf(path, pass_only=False), iter_vcf(tmp, pass_only=False), strict=True):
        if a != b:
            tmp.unlink(missing_ok=True)
            bgzf.gzi_path(tmp).unlink(missing_ok=True)
            raise SystemExit(f"{path.name}: variant {n} differs: {a} != {b}")
        n += 1

    bgzf.gzi_path(tmp).rename(gz.with_name(gz.name + ".gzi"))
    tmp.rename(gz)
    if not keep:
        path.unlink()
    row = {
        "sha256": plain_hash,
        "uncompressed_bytes": before,
        "was_bytes": before,
        "now_bytes": gz.stat().st_size,
        "variants": n,
        "verified": f"{n:,} variants parsed identically and in the same order from both forms",
    }
    print(f"{path.name}: {before / 1e6:.0f} -> {row['now_bytes'] / 1e6:.0f} MB, {n:,} variants identical")
    return row


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", action="append", default=[])
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--skip", action="append", default=[], help="chromosomes to leave alone")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--keep", action="store_true", help="write the blocked file but keep the flat .fa")
    ap.add_argument("--list", action="store_true")
    ap.add_argument(
        "--vcf",
        action="store_true",
        help="compress the split HG002_chr*.vcf files instead of the chromosomes",
    )
    args = ap.parse_args()

    if args.vcf:
        manifest = load_manifest()
        rows = manifest.setdefault("variant_files", {})
        reclaimed = 0
        for p in sorted(REFERENCE.glob("HG002_chr*.vcf")):
            row = compact_vcf(p, keep=args.keep)
            if row is None:
                continue
            reclaimed += row["was_bytes"] - row["now_bytes"]
            rows[p.name] = row
            save_manifest(manifest)
        print(f"reclaimed {reclaimed / 1e9:.2f} GB")
        return 0

    if args.list:
        open_now = held_open()
        for c in chromosomes():
            flat = REFERENCE / f"{c}.fa"
            gz = REFERENCE / f"{c}.fa.gz"
            both = flat.stat().st_size + (gz.stat().st_size if gz.exists() else 0)
            mark = "  HELD OPEN" if {flat.resolve(), gz.resolve()} & open_now else ""
            print(f"{c:14} {both / 1e6:7.0f} MB{mark}")
        return 0

    todo = args.chrom or (chromosomes() if args.all else [])
    todo = [c for c in todo if c not in args.skip]
    if not todo:
        ap.error("name --chrom or --all")

    manifest = load_manifest()
    reclaimed = 0
    for chrom in todo:
        row = compact(chrom, dry_run=args.dry_run, keep=args.keep)
        if row is None:
            continue
        reclaimed += row["was_bytes"] - row["now_bytes"]
        manifest["chromosomes"][chrom] = row
        save_manifest(manifest)
    print(f"reclaimed {reclaimed / 1e9:.2f} GB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
