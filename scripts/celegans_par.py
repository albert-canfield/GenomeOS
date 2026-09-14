# SPDX-License-Identifier: AGPL-3.0-or-later
"""PAR polarity for the first divisions, scored against the published knockout set (area E).

The PAR genes were the only misses the worm program had against Digital Development (Du et al. 2014,
Cell 156:359): four observed founder transformations in par-2 and par-3 that the program could not
produce, because its first two divisions segregated the maternal factors unconditionally. They now
depend on the PAR domain that does the segregating, which is two extra decisions in each founder
layer and no change at all to the wild type:

    div_P0       with PAR-3: the somatic factors (SKN-1, PAL-1) are kept out of AB
    div_P0_no_par3   without it, AB inherits them too and takes the EMS fate
    div_P1       with PAR-2: PIE-1 is retained on the posterior side, in P2
    div_P1_no_par2   without it, neither daughter keeps PIE-1, so P2 takes the EMS fate, presents
                 neither APX-1 nor MOM-2, and ABp and E lose their inductions in turn

This script re-scores both founder layers against the table already distilled in
data/results/celegans_digital_development.json, so it makes no network call, and rewrites that
result with the comparison recomputed from the current programs.

    uv run python scripts/celegans_par.py
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path and not (Path.cwd() / "genomeos").is_dir():
    sys.path.insert(0, str(ROOT))

from genomeos.lang import parse_file  # noqa: E402
from genomeos.organism import digital_development as dd  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = Path("data/results/celegans_digital_development.json")
PROGRAMS = {
    "named_senders": "data/organisms/celegans/embryo.bio",
    "contacts": "data/organisms/celegans/embryo_contacts.bio",
}


def main() -> None:
    t0 = time.time()
    previous = json.loads(RESULT.read_text())
    table = previous["table"]
    runs = {name: dd.compare(parse_file(path), table) for name, path in PROGRAMS.items()}
    main_run = runs["named_senders"]
    out = dict(previous)
    out.update({k: v for k, v in main_run.items() if k != "genes"})
    out["genes"] = main_run["genes"]
    out["table"] = table
    out["programs"] = {
        name: {
            "program": PROGRAMS[name],
            "reproduced": r["reproduced"],
            "missed": r["missed"],
            "not_modelled": r["not_modelled"],
            "rate_over_modelled": r["rate_over_modelled"],
            "missed_changes": [
                {"gene": g["gene"], "cell": o["cell"], "adopts": o["adopts"]}
                for g in r["genes"]
                for o in g["observed"]
                if o["status"] == "missed"
            ],
        }
        for name, r in runs.items()
    }
    out["par_rules"] = {
        "what": "the first two divisions segregate the maternal factors only while the PAR domain "
        "that does the segregating is there",
        "before": "7 of 11 modelled transformations reproduced; the 4 misses were all par-2 and par-3",
        "after": f"{main_run['reproduced']} of {main_run['reproduced'] + main_run['missed']} reproduced",
        "wild_type_unchanged": "PAR-2 and PAR-3 are maternal factors of the zygote, so the conditional "
        "decision is the one that applies in every wild-type run",
        "extra_prediction": "gating the MOM-2 ligand on PIE-1 also predicts E adopting MS in pie-1, "
        "which the published table does not list for that gene; it is a prediction, not a miss",
        "evidence": "Kemphues et al. 1988, Cell 52:311; Bowerman et al. 1993, Cell 74:443; "
        "Boyd et al. 1996, Development 122:3075; Cheeks et al. 2004, Curr Biol 14:851",
    }
    out["seconds"] = round(time.time() - t0, 1)
    print(save_result("celegans_digital_development", out))
    for name, r in runs.items():
        print(name, r["reproduced"], "reproduced,", r["missed"], "missed,", r["not_modelled"], "not modelled")


if __name__ == "__main__":
    main()
