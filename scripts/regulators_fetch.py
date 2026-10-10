# SPDX-License-Identifier: AGPL-3.0-or-later
"""M1b fetch: download the knockdown and control quantification files the registration names.

    uv run python scripts/regulators_fetch.py

Refuses without data/results/regulators_registration.json: the files are fetched only after the test is
registered. Each file is streamed from the ENCODE portal, its md5 checked against the one the
registration pinned from ENCODE's metadata, and kept gzip-compressed and git-ignored under
data/knowledge/regulators/knockdown/. It parses no value: the bytes are hashed and compressed only.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import sys
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos.attribution import regulators as rg  # noqa: E402
from genomeos.results import RESULTS_DIR  # noqa: E402

PORTAL = "https://www.encodeproject.org"
LEDGER = rg.KD_DIR / "fetched.json"


def fetch(f: dict) -> dict:
    out = rg.KD_DIR / f"{f['accession']}.tsv.gz"
    md5, n = hashlib.md5(usedforsecurity=False), 0
    tmp = out.with_suffix(".part")
    req = urllib.request.Request(PORTAL + f["href"])
    with urllib.request.urlopen(req, timeout=300) as r, gzip.open(tmp, "wb", mtime=0) as gz:
        for chunk in iter(lambda: r.read(1 << 20), b""):
            md5.update(chunk)
            n += len(chunk)
            gz.write(chunk)
    if md5.hexdigest() != f["md5sum"]:
        tmp.unlink()
        raise SystemExit(f"refused: {f['accession']} md5 {md5.hexdigest()} is not ENCODE's {f['md5sum']}")
    tmp.rename(out)
    return {"md5": md5.hexdigest(), "bytes": n, "gz_sha256": hashlib.sha256(out.read_bytes()).hexdigest()}


def main() -> int:
    reg_path = RESULTS_DIR / "regulators_registration.json"
    if not reg_path.exists():
        raise SystemExit("refused: register the test before fetching its endpoint")
    reg = json.loads(reg_path.read_text())
    files = {
        x["accession"]: x
        for u in reg["units"]
        for e in u["experiments"]
        for x in e["knockdown_files"] + e["control_files"]
    }
    rg.KD_DIR.mkdir(parents=True, exist_ok=True)
    ledger = json.loads(LEDGER.read_text()) if LEDGER.exists() else {}
    t0 = time.time()
    for i, (acc, f) in enumerate(sorted(files.items()), 1):
        if acc in ledger and (rg.KD_DIR / f"{acc}.tsv.gz").exists():
            continue
        for attempt in range(3):
            try:
                ledger[acc] = fetch(f)
                break
            except SystemExit:
                raise
            except Exception as e:
                if attempt == 2:
                    raise
                print(f"  retry {acc}: {e}", flush=True)
                time.sleep(10)
        LEDGER.write_text(json.dumps(ledger, indent=1, sort_keys=True))
        print(f"{i}/{len(files)} {acc} {ledger[acc]['bytes']} {time.time() - t0:.0f}s", flush=True)
    print("fetched", len(ledger), "files,", sum(v["bytes"] for v in ledger.values()), "bytes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
