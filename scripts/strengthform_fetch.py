# SPDX-License-Identifier: AGPL-3.0-or-later
"""Fetch the one held-out table lane-strengthform reads, and verify its pinned digest.

    uv run python scripts/strengthform_fetch.py

One public GEO supplementary file (GSE240586, Martinez-Ara et al. 2024), about 12.8 MB, into the
git-ignored store data/cache/strengthform. No key, no money, no model request. The digest in
`genomeos.attribution.strengthform.TABLE_SHA256` is what the registration pins and what the run
re-checks, so a file that does not match it is left in place and named rather than used.
"""

from __future__ import annotations

import hashlib
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos.attribution import strengthform as sf  # noqa: E402


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    sf.STORE.mkdir(parents=True, exist_ok=True)
    if not sf.TABLE.exists():
        print("fetching", sf.TABLE_URL)
        with urllib.request.urlopen(sf.TABLE_URL) as response:  # noqa: S310 (a pinned https URL)
            sf.TABLE.write_bytes(response.read())
    got = digest(sf.TABLE)
    print(sf.TABLE, sf.TABLE.stat().st_size, "bytes", got)
    if got != sf.TABLE_SHA256:
        print(f"refused: the digest is {got}, the pinned digest is {sf.TABLE_SHA256}")
        return 1
    print("digest matches the pinned value")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
