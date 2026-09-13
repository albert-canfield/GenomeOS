# SPDX-License-Identifier: AGPL-3.0-or-later
"""Kinetic pathways: curated ODE models from BioModels, run on the in-house SBML engine.

Reactome tells us which reactions a pathway has; it has no rate laws, so `genomeos pathway`
runs it as reachability. BioModels holds hand-curated kinetic models (BIOMD… ids) with rate
laws, parameters and initial conditions for many of the same pathways. This module finds
them by name, keeps the SBML locally (data/knowledge/biomodels, not committed), runs them,
and knocks a species out (initial amount 0, held there) to see what the dynamics do; the
comparison with the untouched run is the kinetic counterpart of Reactome's "reactions lost".

Evidence: the model is curated (BioModels, the paper it reproduces); the run is derived (our
RK4 integration of its equations); a knockout's effect is inferred from that run. Time units
are the model's own.
"""

from __future__ import annotations

import json
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

BIOMODELS = "https://www.ebi.ac.uk/biomodels"
CACHE = Path("data/knowledge/biomodels")
EVIDENCE = "curated: BioModels (manually curated, reproduces the cited paper); derived: in-house RK4 run"


def _get(url: str, accept: str = "application/json", timeout: int = 60) -> bytes:
    req = urllib.request.Request(url, headers={"Accept": accept, "User-Agent": "GenomeOS/0.1 (kinetic)"})
    with urllib.request.urlopen(req, timeout=timeout) as r:  # noqa: S310
        return r.read()


def search(query: str, limit: int = 10, curated_only: bool = True) -> list[dict[str, Any]]:
    """Curated BioModels entries matching a free-text query (a pathway name, a gene, a paper)."""
    url = f"{BIOMODELS}/search?query={urllib.parse.quote(query)}&numResults={max(limit * 3, 10)}&format=json"
    d = json.loads(_get(url))
    out = []
    for m in d.get("models", []):
        mid = m.get("id", "")
        if curated_only and not mid.startswith("BIOMD"):
            continue
        if m.get("format", "SBML") != "SBML":
            continue
        out.append(
            {"id": mid, "name": m.get("name", ""), "submitter": m.get("submitter", ""), "format": "SBML"}
        )
        if len(out) >= limit:
            break
    return out


_STOP = {
    "of",
    "the",
    "by",
    "and",
    "in",
    "to",
    "a",
    "an",
    "via",
    "cascade",
    "pathway",
    "signaling",
    "signalling",
}


def search_pathway(name: str, limit: int = 10) -> tuple[list[dict[str, Any]], str]:
    """A Reactome pathway name as a BioModels query: punctuation out, then ever shorter word sets until
    something curated matches. Returns the hits and the query that found them."""
    import re

    words = [w for w in re.split(r"[^A-Za-z0-9]+", name) if w]
    tries = [" ".join(words)]
    content = [w for w in words if w.lower() not in _STOP and len(w) > 2]
    if content and " ".join(content) not in tries:
        tries.append(" ".join(content))
    for n in (2, 1):
        for i in range(0, max(0, len(content) - n + 1)):
            q = " ".join(content[i : i + n])
            if q not in tries:
                tries.append(q)
    for q in tries:
        try:
            hits = search(q, limit=limit)
        except Exception:  # noqa: BLE001 - a bad query is a 400; try the next
            hits = []
        if hits:
            return hits, q
    return [], name


def fetch(model_id: str, cache_dir: Path = CACHE) -> Path:
    """The model's main SBML file, downloaded once."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"{model_id}.xml"
    if path.exists():
        return path
    info = json.loads(_get(f"{BIOMODELS}/{model_id}?format=json"))
    files = info.get("files", {}).get("main", [])
    if not files:
        raise LookupError(f"{model_id}: no main file on BioModels")
    fname = files[0]["name"]
    path.write_bytes(
        _get(f"{BIOMODELS}/model/download/{model_id}?filename={urllib.parse.quote(fname)}", "*/*")
    )
    (cache_dir / f"{model_id}.json").write_text(
        json.dumps(
            {
                "id": model_id,
                "name": info.get("name", ""),
                "file": fname,
                "publication": info.get("publication", {}),
            }
        )
    )
    return path


def species_matching(model, term: str, accession: str | None = None) -> list[str]:
    """Species whose id or display name contains the term (case-insensitive), or whose MIRIAM annotation
    names the UniProt accession: how a gene symbol lands on a model that calls its species x1, x2, x3."""
    t = term.lower()
    out = []
    for sid in model.species:
        if (
            t in sid.lower()
            or t in model.names.get(sid, "").lower()
            or accession
            and any(u.rstrip("/").endswith(f"/{accession}") for u in model.annotations.get(sid, []))
        ):
            out.append(sid)
    return out


def run(
    path: Path,
    duration: float = 100.0,
    dt: float = 0.01,
    knockout: list[str] | None = None,
    accessions: dict[str, str] | None = None,
) -> dict[str, Any]:
    """Run the model, optionally with species held at zero; per-species final level and peak.
    `accessions` maps a knockout term to a UniProt accession so annotated species match by identity."""
    from genomeos.runtime.sbml import SbmlModel, SbmlRuntime

    model = SbmlModel.from_file(path)
    held = []
    if knockout:
        for term in knockout:
            for sid in species_matching(model, term, (accessions or {}).get(term)):
                model.species[sid] = 0.0
                model.constant_species.add(sid)
                held.append(sid)
    rt = SbmlRuntime(model)
    traj = rt.run(duration=duration, dt=dt, record_every=max(1, int(duration / dt / 400)))
    levels = {}
    for s in traj.species:
        xs = traj.levels[s]
        levels[s] = {
            "name": model.names.get(s, s),
            "initial": round(xs[0], 6),
            "final": round(xs[-1], 6),
            "peak": round(max(xs), 6),
            "peaks": traj.peaks(s),
        }
    return {
        "model": model.id,
        "name": model.name,
        "species": len(model.species),
        "reactions": len(model.reactions),
        "functions": len(model.functions),
        "rate_rules": len(model.rate_rules),
        "duration": duration,
        "dt": dt,
        "knocked_out": held,
        "levels": levels,
        "times": [round(t, 4) for t in traj.times],
        "series": {s: [round(x, 6) for x in traj.levels[s]] for s in traj.species},
        "evidence": EVIDENCE,
    }


def compare(
    baseline: dict[str, Any], knocked: dict[str, Any], min_change: float = 0.1
) -> list[dict[str, Any]]:
    """Species whose final level or peak moved by at least `min_change` of the baseline range: a transient
    that disappears counts as much as a steady state that shifts."""
    out = []
    for s, b in baseline["levels"].items():
        k = knocked["levels"].get(s)
        if not k or s in knocked["knocked_out"]:
            continue
        scale = max(abs(b["peak"]), abs(b["final"]), 1e-9)
        d_final = (k["final"] - b["final"]) / scale
        d_peak = (k["peak"] - b["peak"]) / scale
        delta = d_final if abs(d_final) >= abs(d_peak) else d_peak
        if abs(delta) >= min_change:
            out.append(
                {
                    "species": s,
                    "name": b["name"],
                    "baseline_final": b["final"],
                    "knockout_final": k["final"],
                    "baseline_peak": b["peak"],
                    "knockout_peak": k["peak"],
                    "relative_change": round(delta, 3),
                    "final_change": round(d_final, 3),
                    "peak_change": round(d_peak, 3),
                    "peaks_before": b["peaks"],
                    "peaks_after": k["peaks"],
                }
            )
    out.sort(key=lambda r: -abs(r["relative_change"]))
    return out
