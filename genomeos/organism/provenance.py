# SPDX-License-Identifier: AGPL-3.0-or-later
"""What the worm's results are made from, stated once for their manifests (genomeos/manifest.py).

Item 12 S6 follow-up (lane-contract, 2026-09-28): the C. elegans writers put their results into
data/results without `save_result`, so none carried the provenance contract. They all read the same two
public records through the same two distilled files, so the sources and the common inputs live here and
each writer adds its own parameters, programs, exclusions and partitions.

- The reference lineage: WormWeb's json-celllineage.js (Sulston et al. 1983, Sulston and Horvitz 1977),
  distilled by `genomeos data distil` into data/results/celegans_lineage_cells.json.
- The factor atlas: Ma et al. 2021 (Zenodo 4737593), distilled into
  data/results/celegans_tf_atlas_cells.json (presence per cell) and data/knowledge/celegans/
  atlas_levels.json (levels, built by genomeos.organism.atlas_levels).
"""

from __future__ import annotations

import hashlib
from collections.abc import Iterable
from pathlib import Path
from typing import Any

from genomeos import manifest as mf

from .atlas_levels import LEVELS_FILE
from .reference import KNOWLEDGE as LINEAGE_FILE
from .reference import WORMWEB_URL
from .tf_atlas import CELLS_FILE, FRAME_MIN, PRESENCE_FRACTION, ZENODO_URL

WORMWEB = {
    "accession": "WormWeb C. elegans cell lineage, json-celllineage.js (Sulston et al. 1983; Sulston and "
    "Horvitz 1977)",
    "version": "as fetched and distilled into data/results/celegans_lineage_cells.json (first 0a8aa05, "
    "2026-09-11); the page carries no release number, so the distilled file's sha256 pins it",
    "url": WORMWEB_URL,
    "licence": "CC BY",
}
MA2021 = {
    "accession": "Zenodo 4737593 (Ma et al. 2021, Nat Methods 18:893), Time_series_expression_dataset.zip",
    "version": "Zenodo record 4737593",
    "url": ZENODO_URL,
    "licence": "CC BY 4.0",
}
ASSEMBLY = "n/a: a cell lineage and per-cell factor presence; nothing is placed on a genome"
COORDINATES = "n/a: cells are named by their lineage path, not by a genomic interval"
NO_SPLIT = "n/a: no evaluation split; every cell is read"
#: the worm's hold-out (genomeos.organism.fate_rules.cross_validate and the scripts that follow it)
FOUNDER_HOLDOUT = {
    "held_out_by_founder_sublineage": "rules or thresholds learned on seven of the eight founder "
    "sublineages (fate_rules.SUBLINEAGES: ABal, ABar, ABpl, ABpr, MS, E, C, D) decide the eighth",
    "in_sample": "every score the result marks in sample reads all cells",
}
ORG = Path("data/organisms/celegans")


def bytes_entry(label: str, data: bytes, partition: str | None = None, **extra: Any) -> dict[str, Any]:
    """An input read over the network and not kept (stream, distil, discard): its sha256 and size as
    read, under a label that says where it came from, since no path on disk holds it."""
    return {
        "path": label,
        "sha256": hashlib.sha256(data).hexdigest(),
        "bytes": len(data),
        "partition": partition,
        "kept": False,
        **extra,
    }


def inputs(
    *,
    lineage: bool = True,
    atlas: bool = True,
    levels: bool = False,
    programs: Iterable[str | Path] = (),
    extra: Iterable[str | Path] = (),
) -> list[dict[str, Any]]:
    """The input entries of a worm result, hashed now: call it before a run that rewrites a program, so
    the digest is of what was read."""
    out = []
    if lineage:
        out.append(mf.input_entry(LINEAGE_FILE))
    if atlas:
        out.append(mf.input_entry(CELLS_FILE))
    if levels:
        out.append(mf.input_entry(LEVELS_FILE))
    out += [mf.input_entry(p) for p in programs]
    out += [mf.input_entry(p) for p in extra]
    return out


def manifest(
    entries: list[dict[str, Any]],
    parameters: dict[str, Any],
    *,
    exclusions: list[Any] | None = None,
    partitions: dict[str, Any] | str = NO_SPLIT,
    atlas: bool = True,
) -> dict[str, Any]:
    """A worm result's manifest: the two sources, the inputs hashed by `inputs`, the writer's parameters
    (with the atlas presence rule when the atlas is read), its exclusions and partitions."""
    params = dict(parameters)
    if atlas:
        params.setdefault("atlas_presence_fraction_of_max", PRESENCE_FRACTION)
        params.setdefault("atlas_minutes_per_frame", FRAME_MIN)
    return {
        "sources": [WORMWEB, MA2021] if atlas else [WORMWEB],
        "inputs": entries,
        "assembly": ASSEMBLY,
        "coordinates": COORDINATES,
        "parameters": params,
        "exclusions": list(exclusions or []),
        "partitions": partitions,
    }
