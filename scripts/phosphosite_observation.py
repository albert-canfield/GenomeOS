# SPDX-License-Identifier: AGPL-3.0-or-later
"""Join the curated phosphosites to the Ochoa et al. 2020 reference phosphoproteome.

Stream, distil, discard: the funscoR R data files (about 0.6 MB used) are fetched into a temporary
directory, read, and deleted; only the join to this project's curated sites is written, to
data/results/phosphosite_observation.json. The registration it is scored against is in docs/PROTEIN.md,
"Phosphosite observation (2026-09-27)".

Reading .rda needs the pure-Python ``rdata`` reader, which is not a project dependency:

    uv run --with rdata python scripts/phosphosite_observation.py
"""

from __future__ import annotations

import datetime as dt
import json
import sys
import tempfile
import urllib.request
from pathlib import Path
from statistics import median

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from genomeos.molecules import ptm  # noqa: E402

RAW = "https://raw.githubusercontent.com/evocellnet/funscoR/master/data/feature_spectral_counts.rda"
SOURCE = (
    "Ochoa D, Jarnuczak AF, et al. The functional landscape of the human phosphoproteome. Nat Biotechnol "
    "38:365-373 (2020), doi:10.1038/s41587-019-0344-3; data from the authors' funscoR package "
    "(https://github.com/evocellnet/funscoR, License: LGPL), feature_spectral_counts.rda; "
    "raw data PRIDE PXD012174"
)
CHECK_PRESENT = ("P06748", 125)  # NPM1 S125, constitutive CK2 site
CHECK_ABSENT = ("P04637", 1)  # TP53 M1


def reference() -> dict[tuple[str, int], tuple[str, int, int]]:
    import rdata

    with tempfile.TemporaryDirectory() as tmp:
        f = Path(tmp) / "feature_spectral_counts.rda"
        urllib.request.urlretrieve(RAW, f)
        df = next(iter(rdata.read_rda(str(f)).values()))
    out = {}
    for acc, res, pos, n, psm in zip(
        df["acc"], df["residue"], df["position"], df["Biological_samples"], df["Spectral_Counts"], strict=True
    ):
        out[(str(acc), int(pos))] = (str(res), int(n), int(psm))
    return out


def definitions():
    for p in sorted(ptm.CACHE.glob("*.json")):
        if p.name.startswith("_"):
            continue
        try:
            d = json.loads(p.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        if "sections" in d:
            yield d


def main() -> None:
    ref = reference()
    j = ptm.join_observations(ref, definitions())
    ns = [s[2] for v in j["sites"].values() for s in v]
    present = any(s[0] == CHECK_PRESENT[1] for s in j["sites"].get(CHECK_PRESENT[0], []))
    absent = not any(s[0] == CHECK_ABSENT[1] for s in j["sites"].get(CHECK_ABSENT[0], []))
    ph = j["curated_phospho_sites"]
    matches = j["accession_position_matches"]
    reg_lo, reg_hi = 0.45, 0.70
    frac = j["joined"] / ph if ph else 0.0
    result = {
        "result": "phosphosite_observation",
        "date": dt.date.today().isoformat(),
        "field": "observed_in_cell_types_or_tissues",
        "meaning": ptm.OBSERVATION_MEANING,
        "source": SOURCE,
        "reference_sites": len(ref),
        "curated_sites": j["curated_sites"],
        "curated_phospho_sites": ph,
        "accession_position_matches": matches,
        "residue_disagreements": j["residue_disagreements"],
        "joined": j["joined"],
        "proteins_with_joined_sites": len(j["sites"]),
        "fraction_of_curated_phospho": round(frac, 4),
        "fraction_of_all_curated": round(j["joined"] / j["curated_sites"], 4) if j["curated_sites"] else 0.0,
        "reference_sites_not_joined": len(ref) - j["joined"],
        "other_class_on_reference_position": j["other_class_on_reference_position"],
        "cell_types_or_tissues": {
            "median": median(ns) if ns else None,
            "in_1": sum(1 for n in ns if n == 1),
            "in_10_or_more": sum(1 for n in ns if n >= 10),
            "max": max(ns) if ns else None,
        },
        "registration": {
            "expected_fraction_of_curated_phospho": [reg_lo, reg_hi],
            "within": reg_lo <= frac <= reg_hi,
            "npm1_s125_present": present,
            "tp53_m1_absent": absent,
            "off_acceptor_joined": j["off_acceptor"],
            "residue_disagreement_rate": round(j["residue_disagreements"] / matches, 4) if matches else None,
            "residue_disagreement_under_2pct": (j["residue_disagreements"] / matches < 0.02)
            if matches
            else False,
        },
        "sites_format": (
            "accession -> [[position, residue, observed_in_cell_types_or_tissues, spectral_count], ...]"
        ),
        "sites": j["sites"],
    }
    ptm.OBSERVATION.parent.mkdir(parents=True, exist_ok=True)
    ptm.OBSERVATION.write_text(json.dumps(result, separators=(",", ":")))
    print(json.dumps({k: v for k, v in result.items() if k != "sites"}, indent=1))


if __name__ == "__main__":
    main()
