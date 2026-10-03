#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""N1b's access to the two Replogle pseudobulk files, and nothing else.

Two stages, deliberately apart:

  probe   -- the normalized file's HEADERS and its small obs columns, by HTTP range request. No byte
             of `X` is decoded. This runs BEFORE the registration, because the registration states
             the exact byte ranges the run will ask for and the row counts it will use, and it cannot
             state them without the file's own layout and its own index. What the probe reads is
             design metadata -- which row is non-targeting, which is a core control, how many cells a
             row pooled -- and no outcome.
  rows    -- the exact bytes of the named rows of `X`, one range request per row, each logged with
             its URL, its byte range and its sha256. Row 0 of this file is 32,992 bytes; 585 rows are
             19,300,320 bytes. This runs only AFTER the registration is committed.

Everything written lands under data/cache/, which is never committed: Albert has ruled "data/cache
stays local". A result declares it by path and sha256.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from genomeos.attribution import n1b_calibration as n1b  # noqa: E402

CACHE = ROOT / "data/cache/n1b"
PROBE = CACHE / "probe.json"
RAW_LOCAL = ROOT / "data/cache/n1/K562_gwps_raw_bulk_01.h5ad"
OBS_COLUMNS = ("num_cells_filtered", "core_control", "control_expr")
NON_TARGETING_MARK = "non-targeting"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def probe() -> dict[str, Any]:
    """The normalized file's layout, index and obs columns, by range request. Nothing of `X` decoded."""
    import h5py

    src = n1b.source_file("normalized")
    rf = n1b.RangeFile(src["url"], src["bytes"])
    with h5py.File(io.BufferedReader(rf, buffer_size=n1b.BLOCK), "r") as h:
        layout = n1b.x_layout(h)
        index = [x.decode() if isinstance(x, bytes) else str(x) for x in h["obs/gene_transcript"][()]]
        var_ids = [x.decode() if isinstance(x, bytes) else str(x) for x in h["var/gene_id"][()]]
        obs = {c: [float(v) for v in h[f"obs/{c}"][()]] for c in OBS_COLUMNS}
    nt = [i for i, s in enumerate(index) if NON_TARGETING_MARK in s]
    usable = n1b.usable_rows([obs["num_cells_filtered"][i] for i in nt])
    rec = {
        "stage": "probe",
        "what_was_read": "the normalized file's superblock and object headers, obs/gene_transcript, "
        "var/gene_id, and three obs columns. No byte of X was decoded.",
        "source": {**src, "md5_whole_file": n1b.MD5_NOT_VERIFIABLE},
        "x_layout": layout,
        "rows_total": len(index),
        "genes_total": len(var_ids),
        "non_targeting_rows": nt,
        "non_targeting_count": len(nt),
        "core_control_count_overall": int(sum(1 for v in obs["core_control"] if v)),
        "core_control_count_among_non_targeting": int(sum(1 for i in nt if obs["core_control"][i])),
        "usable_non_targeting": len(usable["rows_used"]),
        "unusable_non_targeting": len(usable["rows_dropped"]),
        "control_expr_finite_among_non_targeting": int(
            sum(1 for i in nt if n1b._finite(obs["control_expr"][i]))  # noqa: SLF001
        ),
        "obs": {c: [obs[c][i] for i in nt] for c in OBS_COLUMNS},
        "bytes_fetched": rf.fetched,
        "fetch_log": rf.log,
    }
    CACHE.mkdir(parents=True, exist_ok=True)
    PROBE.write_text(json.dumps(rec, indent=1, sort_keys=True))
    return rec


def row_ranges(layout: dict[str, Any], rows: list[int]) -> dict[str, Any]:
    """The exact byte range of each named row, stated so a registration can carry it before the read."""
    base, stride = layout["offset"], layout["row_stride_bytes"]
    return {
        "dataset": "X",
        "base_offset": base,
        "row_stride_bytes": stride,
        "rows": len(rows),
        "first_range": f"bytes={base + rows[0] * stride}-{base + (rows[0] + 1) * stride - 1}",
        "last_range": f"bytes={base + rows[-1] * stride}-{base + (rows[-1] + 1) * stride - 1}",
        "total_bytes": len(rows) * stride,
        "requests": len(rows),
    }


def fetch_rows(rows: list[int], layout: dict[str, Any]) -> tuple[Any, list[dict[str, Any]]]:
    """The named rows of normalized `X`, by one exact range request each. Returns (array, fetch log)."""
    import numpy as np

    src = n1b.source_file("normalized")
    rf = n1b.RangeFile(src["url"], src["bytes"])
    cols = layout["cols"]
    out = np.empty((len(rows), cols), dtype=np.float32)
    base, stride = layout["offset"], layout["row_stride_bytes"]
    for p, r in enumerate(rows):
        blob = rf._get(base + r * stride, base + (r + 1) * stride - 1, f"X row {r}")  # noqa: SLF001
        out[p] = np.frombuffer(blob, dtype="<f4", count=cols)
        del blob
    return out, rf.log


def raw_expression(rows: list[int], cell_counts: list[float], block: int = 64) -> Any:
    """Per gene, the cell-weighted mean of RAW X over the given rows, read from the local raw file in
    blocks of `block` rows so no whole matrix is opened. This is the stratifier's measured value."""
    import h5py
    import numpy as np

    total = float(sum(cell_counts))
    acc = None
    with h5py.File(RAW_LOCAL, "r") as h:
        ds = h["X"]
        for s in range(0, len(rows), block):
            idx = rows[s : s + block]
            chunk = np.asarray(ds[idx], dtype=np.float64)
            w = np.asarray(cell_counts[s : s + block], dtype=np.float64)
            part = (chunk * w[:, None]).sum(axis=0)
            acc = part if acc is None else acc + part
            del chunk, part
    return acc / total


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("stage", choices=["probe", "ranges"])
    args = ap.parse_args()
    if args.stage == "probe":
        rec = probe()
        print(
            f"probe: {rec['rows_total']} rows, {rec['genes_total']} genes, "
            f"{rec['non_targeting_count']} non-targeting, "
            f"{rec['core_control_count_among_non_targeting']} core, "
            f"{rec['usable_non_targeting']} usable, "
            f"control_expr finite on {rec['control_expr_finite_among_non_targeting']} of them; "
            f"{rec['bytes_fetched']} bytes in {len(rec['fetch_log'])} range requests"
        )
        print(f"X at offset {rec['x_layout']['offset']}, stride {rec['x_layout']['row_stride_bytes']}")
        return 0
    rec = json.loads(PROBE.read_text())
    print(json.dumps(row_ranges(rec["x_layout"], rec["non_targeting_rows"]), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
