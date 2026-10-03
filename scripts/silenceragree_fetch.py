#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Fetch the ReSE screen-validated silencer records from NCBI Gene, free and without a key.

Writes data/cache/silenceragree/rese_records.json, which is NEVER committed (data/cache is
git-ignored). The registration and the run declare it by path and sha256.

NCBI eutils takes no key for this rate and charges nothing. Zero AlphaGenome requests, no money.
"""

from __future__ import annotations

import hashlib
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
#: The one query. Every record NCBI returns for it is written; the selection to the joinable set is
#: the module's job, not the fetcher's, so the cache keeps what was rejected and why is auditable.
QUERY = "ReSE[All Fields] AND human[orgn]"
OUT = ROOT / "data/cache/silenceragree/rese_records.json"
PAGE = 400


def _get(url: str, tries: int = 5) -> str:
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as resp:  # noqa: S310
                return resp.read().decode()
        except Exception as exc:  # noqa: BLE001
            if attempt == tries - 1:
                raise
            print(f"retry {attempt}: {exc}", file=sys.stderr)
            time.sleep(3)
    raise AssertionError("unreachable")


def main() -> int:
    term = urllib.parse.quote(QUERY)
    head = _get(f"{BASE}esearch.fcgi?db=gene&term={term}&usehistory=y&retmax=0")
    env = re.search(r"<WebEnv>(.*?)</WebEnv>", head)
    key = re.search(r"<QueryKey>(.*?)</QueryKey>", head)
    total = re.search(r"<Count>(\d+)</Count>", head)
    if not (env and key and total):
        raise SystemExit(f"esearch did not answer with a history handle: {head[:400]}")
    n = int(total.group(1))
    print(f"query {QUERY!r} -> {n} records")
    records: dict[str, dict] = {}
    for start in range(0, n, PAGE):
        url = (
            f"{BASE}esummary.fcgi?db=gene&query_key={key.group(1)}&WebEnv={env.group(1)}"
            f"&retstart={start}&retmax={PAGE}&retmode=json"
        )
        res = json.loads(_get(url))["result"]
        for uid in res["uids"]:
            records[uid] = res[uid]
        print(f"  fetched {len(records)}/{n}")
        time.sleep(0.4)
    if len(records) != n:
        raise SystemExit(f"esearch said {n} records and esummary gave {len(records)}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({"query": QUERY, "count": n, "records": records}, sort_keys=True))
    print(f"wrote {OUT.relative_to(ROOT)}")
    print(f"sha256 {hashlib.sha256(OUT.read_bytes()).hexdigest()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
