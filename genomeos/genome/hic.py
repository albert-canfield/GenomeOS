# SPDX-License-Identifier: AGPL-3.0-or-later
"""Measured Hi-C boundaries from the 4D Nucleome portal, held against the CTCF-only nodes.

4DN publishes boundary calls (insulation-based, BED) for hundreds of released GRCh38 Hi-C and Micro-C
experiments; the search is open, the downloads need a free account key (FOURDN_KEY and FOURDN_SECRET in
.env, like the AlphaGenome key). This module finds the boundary files of a biosource, streams them once
into per-chromosome files under data/knowledge/hic (local), and compares the measured boundaries with
the CTCF-only ones the way the predicted contact map was compared, random control included. Measured
evidence at last for the node edges; without a key the feature loads disabled and says so.
"""

from __future__ import annotations

import base64
import gzip
import io
import json
import os
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

FOURDN = "https://data.4dnucleome.org"
KNOWLEDGE = Path("data/knowledge/hic")
KEY_VAR, SECRET_VAR = "FOURDN_KEY", "FOURDN_SECRET"
EVIDENCE = "experimental: 4D Nucleome boundary calls (insulation score) on released GRCh38 Hi-C and Micro-C"
HOW = (
    "create a free 4DN account, generate an access key at https://data.4dnucleome.org and put "
    "FOURDN_KEY=... and FOURDN_SECRET=... in the git-ignored .env"
)


def credentials() -> tuple[str, str] | None:
    from genomeos.predict.alphagenome_adapter import _dotenv_key

    key = os.environ.get(KEY_VAR) or _dotenv_key(KEY_VAR)
    secret = os.environ.get(SECRET_VAR) or _dotenv_key(SECRET_VAR)
    return (key, secret) if key and secret else None


def status() -> dict[str, Any]:
    return {
        "name": "4DN Hi-C boundaries",
        "enabled": credentials() is not None,
        "how": HOW,
        "evidence": EVIDENCE,
    }


class _DropAuthOnRedirect(urllib.request.HTTPRedirectHandler):
    """4DN answers an authenticated download with a redirect to S3; S3 refuses a request that still carries
    the portal's Authorization header, so the header is dropped when the redirect leaves the portal's host."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):  # noqa: PLR0913
        new = super().redirect_request(req, fp, code, msg, headers, newurl)
        if (
            new is not None
            and urllib.parse.urlparse(newurl).netloc != urllib.parse.urlparse(req.full_url).netloc
        ):
            new.remove_header("Authorization")
        return new


def _get(
    url: str, auth: tuple[str, str] | None = None, accept: str = "application/json", timeout: int = 120
) -> bytes:
    headers = {"Accept": accept, "User-Agent": "GenomeOS/0.1 (stream)"}
    if auth:
        token = base64.b64encode(f"{auth[0]}:{auth[1]}".encode()).decode()
        headers["Authorization"] = f"Basic {token}"
    req = urllib.request.Request(url, headers=headers)
    opener = urllib.request.build_opener(_DropAuthOnRedirect())
    with opener.open(req, timeout=timeout) as r:  # noqa: S310
        return r.read()


def find_boundary_files(biosource: str, limit: int = 400) -> list[dict[str, Any]]:
    """Released GRCh38 boundary files of a biosource (the search needs no key)."""
    url = (
        f"{FOURDN}/search/?type=FileProcessed&genome_assembly=GRCh38&status=released&file_type=boundaries"
        f"&limit={limit}&format=json"
    )
    d = json.loads(_get(url))
    out = []
    for f in d.get("@graph", []):
        info = f.get("track_and_facet_info") or {}
        if (info.get("biosource_name") or "") != biosource:
            continue
        out.append(
            {
                "accession": f.get("accession"),
                "href": FOURDN + f["href"],
                "size": f.get("file_size"),
                "experiment_type": info.get("experiment_type"),
                "assay": info.get("assay_info"),
                "biosource": biosource,
            }
        )
    out.sort(key=lambda x: -(x["size"] or 0))  # deepest call set first
    return out


def fetch_boundaries(biosource: str, knowledge: Path = KNOWLEDGE, progress=None) -> dict[str, Any]:
    """Stream the biosource's largest boundary file once (needs the key) into per-chromosome files."""
    slug = biosource.replace(" ", "_").replace("/", "_")
    d = knowledge / slug
    manifest = d / "manifest.json"
    if manifest.exists():
        return json.loads(manifest.read_text())
    auth = credentials()
    if auth is None:
        raise PermissionError(f"4DN downloads need an account key: {HOW}")
    files = find_boundary_files(biosource)
    if not files:
        raise LookupError(f"no released GRCh38 boundary file on 4DN for biosource {biosource!r}")
    chosen = files[0]
    raw = _get(chosen["href"], auth, accept="*/*", timeout=600)
    text = gzip.decompress(raw).decode() if raw[:2] == b"\x1f\x8b" else raw.decode()
    d.mkdir(parents=True, exist_ok=True)
    counts: dict[str, int] = {}
    handles: dict[str, Any] = {}
    try:
        for line in io.StringIO(text):
            if line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            chrom = f[0] if f[0].startswith("chr") else f"chr{f[0]}"
            if chrom not in handles:
                handles[chrom] = open(d / f"{chrom}.bed", "w")  # noqa: SIM115
            handles[chrom].write(f"{chrom}\t{f[1]}\t{f[2]}\n")
            counts[chrom] = counts.get(chrom, 0) + 1
    finally:
        for h in handles.values():
            h.close()
    out = {**chosen, "chromosomes": counts, "evidence": EVIDENCE}
    manifest.write_text(json.dumps(out, indent=1))
    if progress:
        progress(
            f"{biosource}: {sum(counts.values()):,} boundaries on {len(counts)} chromosomes "
            f"({chosen['accession']})"
        )
    return out


def load_boundaries(biosource: str, chrom: str, knowledge: Path = KNOWLEDGE) -> list[int]:
    """Measured boundary positions (interval centres) on a chromosome, sorted."""
    slug = biosource.replace(" ", "_").replace("/", "_")
    p = knowledge / slug / f"{chrom}.bed"
    if not p.exists():
        return []
    out = []
    with p.open() as fh:
        for line in fh:
            f = line.split("\t")
            out.append((int(f[1]) + int(f[2])) // 2)
    return sorted(out)


def compare_with_inferred(inferred: list[int], measured: list[int], length: int) -> dict[str, Any]:
    from genomeos.predict.contact_maps import compare, random_control

    c = compare(inferred, measured)
    c["random_control"] = random_control(inferred, measured, length)
    return c
