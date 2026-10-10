# SPDX-License-Identifier: AGPL-3.0-or-later
"""Area E's fourth number in a program with no lookup: how many of the worm's 555 embryonic terminal
fates are decided, and decided right, when none of the observed lineage's `fate_<cell>` decisions is in
the program at all.

WHY THIS AND NOT A BETTER READ. Area E publishes four numbers together: 522 with the lookup as fallback,
381 honest, 136 free, **0 decided without the lookup having written the fate first** (d5215bf). The 0 is
structural: every factor rule in fates.bio is guarded on a terminal `cell_type`, and only the lookup ever
sets one, so a rule can confirm a fate and never decide one. The terminality gate (de7a135, 9d266a3)
then showed that the factors cannot say *when* a cell is terminal well enough to lift the guard. So the
next number has to come from somewhere other than the terminal factor rules: here, from the founder
decisions, which since 9749b28 read the signals at contacts rather than at named cells, through the three
founders whose whole clone is one tissue (founder_tissues.bio).

WHAT THE PROGRAM IS. data/organisms/celegans/embryo_nolookup.bio: embryo_contacts.bio with the lineage's
555 terminal-fate decisions removed (lineage_embryo_tree.bio, derived by `--write-tree`), the three
founder-tissue rules added, fates.bio left out, and the contact table replaced by the MEASURED one
(contacts_embryo_cshaper.tsv: Cao et al. 2020, CShaper, Supplementary Data 5, the 479 contacts
reproducible in all 17 embryos; scripts/celegans_cshaper_contacts.py). No runtime change was needed:
contact-at-birth reading (`mode: contact; reads: amount`, `Body.neighbours`, the re-resolve of a
newborn's neighbours in `Body._add`) is BioLang v0.4 §7.4 and was already in genomeos/runtime/body.py.

PRE-REGISTRATION. `REGISTRATION` below is written to data/results/celegans_nolookup_registration.json by
`--register` and committed before the program is run for any fate; `run` refuses unless that file is
committed, equal to `REGISTRATION`, and the program's files hash as registered. Before registering, the
program was grown once to check it runs (1,439 cells and 110 deaths by 800 min, no unknowns) and its own
count and death asserts were run by `bio test`; the shipped embryo_contacts.bio was grown with the
timers at their means and with seeds 0 and 7 to check that V2 can hold (1,439 of 1,439 reference cells
matched each time). No fate, founder type or tissue count was read from any of those runs.

The contact partners of the founder cells WERE read from the CShaper sheets before this was written
(ABa touches ABp and EMS only, in every sheet; ABp, EMS, ABpl and E touch P2 in Sheet4), so the
founder-identity predictions are expected rather than blind. The counts, the seeds, the controls and the
knockouts' terminal tissues were not run; their predictions are derived from the rules as written, and
the PIE-1 one was corrected once while writing, before any run, when the derivation was redone and
showed that `presents_mom2` needs PIE-1 in P2.

    uv run python scripts/celegans_nolookup.py --write-tree   # lineage_embryo_tree.bio
    uv run python scripts/celegans_nolookup.py --register      # the registration, then commit it
    uv run python scripts/celegans_nolookup.py                 # the run

No model calls, no network, under a minute on one core.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import tempfile
from collections import Counter
from datetime import date
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path and not (Path.cwd() / "genomeos").is_dir():
    sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.lang import parse_file  # noqa: E402
from genomeos.organism import fate_rules as fr  # noqa: E402
from genomeos.organism import provenance as pv  # noqa: E402
from genomeos.organism.diff import compare  # noqa: E402
from genomeos.organism.reference import ReferenceLineage  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402
from genomeos.runtime.body import Body  # noqa: E402

ORG = Path("data/organisms/celegans")
PROGRAM = ORG / "embryo_nolookup.bio"
LINEAGE = ORG / "lineage_embryo.bio"
TREE = ORG / "lineage_embryo_tree.bio"
MEASURED = ORG / "contacts_embryo_cshaper.tsv"
DETECTED = ORG / "contacts_embryo_cshaper_detected.tsv"
CURATED = ORG / "contacts_embryo.tsv"
REGISTRATION_NAME = "celegans_nolookup_registration"
RESULT_NAME = "celegans_nolookup"
HATCH = 800.0
FOUNDER_UNTIL = 250.0
SEEDS = list(range(20))
OWN_CODE = (
    "scripts/celegans_nolookup.py",
    "scripts/celegans_cshaper_contacts.py",
    "tests/test_celegans_nolookup.py",
)
#: the files the run must find exactly as registered
FROZEN = (
    str(PROGRAM),
    str(TREE),
    str(ORG / "founder_tissues.bio"),
    str(ORG / "founders_contacts.bio"),
    str(MEASURED),
    str(DETECTED),
    str(CURATED),
    str(LINEAGE),
    str(ORG / "lineage_larva.bio"),
    str(ORG / "timers.bio"),
    str(ORG / "reader.bio"),
    str(ORG / "cell_types.bio"),
)
TISSUE_RULES = ("tissue_E", "tissue_D", "tissue_germline")
FOUNDER_DECISIONS = (
    "fate_ABp",
    "fate_ABa",
    "fate_EMS_polarised",
    "fate_EMS",
    "fate_MS",
    "fate_E",
    "fate_C",
    "fate_D",
    "germline",
)
CSHAPER = {
    "accession": "Cao et al. 2020, Nat Commun 11:6254, doi:10.1038/s41467-020-19863-x (CShaper), "
    "Supplementary Data 5 (41467_2020_19863_MOESM6_ESM.xls)",
    "version": "sha256 519c34349af8e8f6cbd675c665f9a880c99fe1ad482b469067f4a6078a855c99, 461,824 bytes, "
    "distilled by scripts/celegans_cshaper_contacts.py",
    "url": "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-19863-x/MediaObjects/"
    "41467_2020_19863_MOESM6_ESM.xls",
    "licence": "CC BY 4.0",
}

# ---- the registration -----------------------------------------------------------------------------

WT_COUNTS = {"Intestine": 34, "Muscle": 20, "Reproductive": 2, "decided": 56}
REGISTRATION: dict[str, Any] = {
    "question": "Of the 555 embryonic terminal cells of the reference lineage, how many end with a "
    "terminal type, and the right one, in a program that contains none of the lookup's 555 fate_<cell> "
    "decisions?",
    "program": str(PROGRAM),
    "what_moves_off_the_lookup": {
        "tissue_E": {
            "cells": "E's 34 terminal descendants, all Intestine in the reference",
            "founder_decided_by": "SKN-1 present and POP-1 absent (fate_E); POP-1 is absent in E only "
            "because EMS was polarised by MOM-2 read at the measured EMS-P2 contact (fate_EMS_polarised, "
            "div_EMS_polarised)",
            "reads_a_measured_contact": True,
            "circular": "yes, as a wild-type prediction: 'E's clone is intestine' is Sulston 1983's "
            "observation, the record the 555 are scored against",
        },
        "tissue_D": {
            "cells": "D's 20 terminal descendants, all Muscle in the reference",
            "founder_decided_by": "the name D with PAL-1 present and SKN-1 absent, both segregated by "
            "named divisions (fate_D)",
            "reads_a_measured_contact": False,
            "circular": "yes, as a wild-type prediction: 'D's clone is body-wall muscle' is Sulston 1983's "
            "observation",
        },
        "tissue_germline": {
            "cells": "P4's two daughters, Z3 (P4a) and Z2 (P4p), Reproductive in the reference",
            "founder_decided_by": "the name P4 with PIE-1 present (germline)",
            "reads_a_measured_contact": False,
            "circular": "yes, and it is also a translation of the program's GermCell into the "
            "reference's word for the germ line",
        },
        "not_moved": "every AB, MS and C terminal cell: none of those founders' clones is one tissue, "
        "and a rule naming a sub-founder would be the lookup under another name",
    },
    "definitions": {
        "denominator": "the 555 embryonic terminal cells: reference cells that are terminal, do not die, "
        "and are born before 800 min (fate_rules.embryonic_terminal)",
        "decided_without_the_lookup": "a cell of the 555 whose cell_type at 800 min is one of the 17 "
        "terminal types, in a program in which no decision is fate_<cell> for any of the 555 (asserted "
        "before counting; the run is void otherwise)",
        "right": "that type equals the reference type; wrong otherwise",
        "attribution": "walking from the cell up its parents in the run, the first cell on which a "
        "differentiate decision to that type fired; every decided cell must be attributed to one of "
        "tissue_E, tissue_D, tissue_germline, or the run is void",
        "primary_arm": "embryo_nolookup.bio as committed, every timer at its mean, grown to 800 min",
    },
    "predictions": {
        "P1 primary": "decided 56 (tissue_E 34, tissue_D 20, tissue_germline 2), right 56, wrong 0",
        "P2 seeds": "the same 56 / 56 / 0 with the timers drawn, seeds 0 to 19, every seed",
        "P3 founders": "the eight founder types under the measured table are the expected ones: ABa "
        "ABaPrecursor, ABp ABpPrecursor, EMS EMSPolarised, MS MSPrecursor, E EPrecursor, C CPrecursor, "
        "D DPrecursor, P4 GermCell (8 of 8)",
        "P4 other tables": "the curated table (contacts_embryo.tsv) and the loosest measured reading "
        "(Sheet1, 3,011 detected pairs) give the same eight founder types and the same 56 / 56 / 0",
        "P5 no table": "with no contact table nothing is induced: ABp takes ABaPrecursor, E takes "
        "MSPrecursor, Intestine 0, decided 22 (tissue_D 20, tissue_germline 2), right 22",
        "P6 rearranged": "ABa and ABp swapped in the measured table alone: ABa takes ABpPrecursor and "
        "ABp takes ABaPrecursor",
        "P7 knockouts, founder level": "the eight knockouts of scripts/celegans_contacts.py give their "
        "published founder transformations under the measured table: 8 of 8",
        "P8 knockouts, terminal level": "Intestine / Muscle / Reproductive / decided among the 555 are "
        "POP-1 111/20/2/133, SKN-1 0/20/2/22, PIE-1 0/20/0/20, APX-1 and GLP-1 34/20/2/56 (as wild "
        "type), PAL-1 34/0/2/36, MOM-2 and MOM-5 0/20/2/22",
    },
    "published_terminal_phenotypes": {
        "POP-1": "MS adopts the E fate: extra intestine, no MS-derived pharynx or muscle (Lin, Thorpe & "
        "Priess 1995, Cell 83:599). Registered reading: reproduced in direction if MS's terminal cells "
        "are Intestine; the count (111) is NOT comparable, because MS's divisions are still read from "
        "the lineage by name, where a pop-1 MS makes an E-like lineage",
        "SKN-1": "no intestine and no pharynx; EMS daughters make C-like hypodermis and muscle (Bowerman "
        "et al. 1992, Cell 68:1061). Registered reading: reproduced if Intestine is 0; the extra C-like "
        "tissue is not decidable here (no C tissue rule)",
        "PIE-1": "P2 takes the EMS fate: extra pharynx and intestine and no germ line (Mello et al. 1992, "
        "Cell 70:163). Registered reading: 'no germ line' reproduced if Reproductive is 0; 'extra "
        "intestine' NOT reproduced, predicted in advance and in the opposite direction: the program "
        "loses E's intestine too, because presents_mom2 requires PIE-1 in P2, so without PIE-1 no Wnt "
        "reaches EMS at the contact and E takes the MS fate; and the MS/E founder rules are guarded on "
        "the names MS and E, so a P2 turned EMS-like still divides into cells named C and P3, and D "
        "keeps its muscle",
        "APX-1": "ABp becomes ABa-like (Mello et al. 1994, Cell 77:95). Registered reading: a founder-"
        "level phenotype only; no AB fate is decided without the lookup, so nothing terminal can move",
        "GLP-1": "ABp becomes ABa-like (Priess, Schnabel & Schnabel 1987, Cell 51:601). As APX-1",
        "PAL-1": "C and D lose their fates: no D muscle (Hunter & Kenyon 1996, Cell 87:217). Registered "
        "reading: reproduced if Muscle is 0; C's tissues are not decided either way",
        "MOM-2": "E takes the MS fate, no intestine (Thorpe et al. 1997, Cell 90:695; Rocheleau et al. "
        "1997, Cell 90:707). Registered reading: reproduced if Intestine is 0",
        "MOM-5": "E takes the MS fate (Rocheleau et al. 1997). The published RNAi phenotype is partially "
        "penetrant; the program's knockout is complete, so only the direction is compared",
    },
    "what_is_circular": [
        "every one of the expected 56 wild-type fates: each founder-tissue rule was written from the "
        "clonal observation in the same Sulston lineage the fates are scored against",
        "every founder-level knockout phenotype: the founder rules cite these same papers",
        "the terminal-level knockout outcomes follow from those two circular rule sets; a pass is "
        "consistency, not prediction, and a fail is a limit of the rules",
        "what is not circular: the contact topology itself (CShaper, a segmentation study made for "
        "another purpose), and so whether the measured contacts deliver the founder signals to the cells "
        "the genetics says they reach",
    ],
    "void_if": [
        "V1 any of the 555 has a fate_<cell> decision in the program",
        "V2 the primary arm does not grow the reference tree to 800 min: the reference cells born by "
        "800 min must all be matched by name with no parent mismatch",
        "V3 genomeos.organism.diff.compare counts a different number of correct terminal fates than "
        "this instrument's right count on the primary arm (two instruments disagree)",
        "V4 a decided cell cannot be attributed to tissue_E, tissue_D or tissue_germline",
    ],
    "falsifiers": {
        "measured contacts": "any of the eight founder types differs under the measured table (P3): "
        "the measured geometry does not deliver the signals as the program writes them",
        "contact dependence": "P5 fails with E still Intestine: the 34 do not depend on a contact and "
        "that claim is withdrawn",
        "count": "P1 off: the founder rules do not carry the clones as written",
    },
    "wording": {
        "P1 holds": "56 of 555 embryonic terminal fates are decided without the lookup, 56 right and 0 "
        "wrong, in embryo_nolookup.bio; all 56 are circular as wild-type predictions (each founder-"
        "tissue rule was written from Sulston's clonal observation, the record they are scored against); "
        "34 follow a Wnt read at a measured contact (Cao et al. 2020), 22 follow maternal factors "
        "segregated by named divisions. The shipped embryo_factors.bio is unchanged and its count "
        "stays 0",
        "P1 fails": "N of 555 decided without the lookup, R right, W wrong; stated with which clone "
        "departs from the registered 34 / 20 / 2",
        "seeds vary": "the count is reported as its range over seeds 0 to 19, and the primary number as "
        "the means run only",
    },
}

# ---- the founder knockouts, as scripts/celegans_contacts.py states them ------------------------------

KNOCKOUT_FOUNDERS: dict[str, dict[str, str]] = {
    "APX-1": {"ABa": "ABaPrecursor", "ABp": "ABaPrecursor"},
    "GLP-1": {"ABa": "ABaPrecursor", "ABp": "ABaPrecursor"},
    "MOM-2": {"MS": "MSPrecursor", "E": "MSPrecursor"},
    "MOM-5": {"MS": "MSPrecursor", "E": "MSPrecursor"},
    "PIE-1": {"P2": "EMSPrecursor", "ABp": "ABaPrecursor"},
    "SKN-1": {"MS": "CPrecursor", "E": "CPrecursor"},
    "POP-1": {"MS": "EPrecursor", "E": "EPrecursor"},
    "PAL-1": {"C": "Blastomere", "D": "Blastomere"},
}
EXPECTED_FOUNDERS = {
    "ABa": "ABaPrecursor",
    "ABp": "ABpPrecursor",
    "EMS": "EMSPolarised",
    "MS": "MSPrecursor",
    "E": "EPrecursor",
    "C": "CPrecursor",
    "D": "DPrecursor",
    "P4": "GermCell",
}
KNOCKOUT_TERMINAL = {
    "POP-1": {"Intestine": 111, "Muscle": 20, "Reproductive": 2, "decided": 133},
    "SKN-1": {"Intestine": 0, "Muscle": 20, "Reproductive": 2, "decided": 22},
    "PIE-1": {"Intestine": 0, "Muscle": 20, "Reproductive": 0, "decided": 20},
    "APX-1": dict(WT_COUNTS),
    "GLP-1": dict(WT_COUNTS),
    "PAL-1": {"Intestine": 34, "Muscle": 0, "Reproductive": 2, "decided": 36},
    "MOM-2": {"Intestine": 0, "Muscle": 20, "Reproductive": 2, "decided": 22},
    "MOM-5": {"Intestine": 0, "Muscle": 20, "Reproductive": 2, "decided": 22},
}
#: the published phenotype at the terminal level, as the registered reading of each one says to test it
PUBLISHED_TEST = {
    "POP-1": lambda c, by: by.get("MS", {}).get("Intestine", 0) > 0,
    "SKN-1": lambda c, by: c["Intestine"] == 0,
    "PIE-1": lambda c, by: c["Reproductive"] == 0 and c["Intestine"] > WT_COUNTS["Intestine"],
    "PAL-1": lambda c, by: by.get("D", {}).get("Muscle", 0) == 0,
    "MOM-2": lambda c, by: c["Intestine"] == 0,
    "MOM-5": lambda c, by: c["Intestine"] == 0,
}
REGISTRATION["predictions"]["P8 detail"] = KNOCKOUT_TERMINAL
REGISTRATION["predictions"]["P8 published terminal phenotype reproduced"] = {
    "POP-1": True,
    "SKN-1": True,
    "PIE-1": "False: 'no germ line' holds, 'extra intestine' does not (Intestine 0)",
    "APX-1": "not testable at the terminal level",
    "GLP-1": "not testable at the terminal level",
    "PAL-1": True,
    "MOM-2": True,
    "MOM-5": "True in direction; penetrance not compared",
}

# ---- the derived lineage ----------------------------------------------------------------------------

TREE_HEADER = (
    "# The observed embryonic lineage WITHOUT its 555 terminal fates (area E): lineage_embryo.bio",
    "# with every `decision fate_<cell>` line removed and nothing else changed -- the divisions, the",
    "# daughters' names and the programmed deaths are Sulston's, as before, and no cell is told what it",
    "# becomes. Derived by `uv run python scripts/celegans_nolookup.py --write-tree`;",
    "# tests/test_celegans_nolookup.py checks the file is exactly that derivation of the lineage it",
    "# comes from, so the two cannot drift apart.",
)


def tree_text(lineage: str) -> str:
    """lineage_embryo.bio without its terminal-fate decisions, under its own module name."""
    out = list(TREE_HEADER)
    for line in lineage.splitlines():
        if line.startswith("decision fate_"):
            continue
        if line.startswith("module organism.celegans.lineage"):
            line = "module organism.celegans.lineage_tree"
        out.append(line)
    return "\n".join(out) + "\n"


# ---- the instrument ---------------------------------------------------------------------------------


def terminal_cells(ref: ReferenceLineage) -> dict[str, str]:
    """The 555, with their reference type."""
    return {c.id: c.cell_type for c in fr.embryonic_terminal(ref)}


def founder_of(ref: ReferenceLineage, cid: str) -> str:
    """The founder a reference cell descends from (AB, MS, E, C, D, P4), by its parents."""
    x = ref.cells.get(cid)
    while x is not None:
        if x.id in ("AB", "MS", "E", "C", "D", "P4"):
            return x.id
        x = ref.cells.get(x.parent) if x.parent else None
    return "none"


def lookup_fates_present(module, cells: dict[str, str]) -> list[str]:  # noqa: ANN001
    """V1: the ids of any fate_<cell> decision for a cell of the 555."""
    names = {f"fate_{c}" for c in cells}
    return sorted(d.id for d in module.decisions if d.id in names)


def attribute(body: Body, cid: str, targets: dict[str, str]) -> str | None:
    """The differentiate decision that set this cell's type: the first cell up its parents on which a
    decision to that type fired."""
    c = body.cells.get(cid)
    want = c.cell_type if c is not None else None
    while c is not None:
        for did in reversed(c.fired):
            if targets.get(did) == want:
                return did
        c = body.cells.get(c.parent) if c.parent else None
    return None


def decided(body: Body, cells: dict[str, str], ref: ReferenceLineage, targets: dict[str, str]) -> dict:
    """The count, cell by cell: decided, right, wrong, by rule and by founder clone."""
    terminal_types = set(fr.TERMINAL_TYPES.split("|"))
    out: dict[str, Any] = {"cells": len(cells), "decided": 0, "right": 0, "wrong": 0}
    by_rule: Counter[str] = Counter()
    by_type: Counter[str] = Counter()
    by_founder: dict[str, Counter[str]] = {}
    wrong: Counter[str] = Counter()
    unattributed: list[str] = []
    for cid, want in sorted(cells.items()):
        c = body.cells.get(cid)
        got = c.cell_type if c is not None else None
        if got not in terminal_types:
            continue
        out["decided"] += 1
        by_type[got] += 1
        by_founder.setdefault(founder_of(ref, cid), Counter())[got] += 1
        rule = attribute(body, cid, targets)
        if rule is None:
            unattributed.append(cid)
        by_rule[rule or "none"] += 1
        if got == want:
            out["right"] += 1
        else:
            out["wrong"] += 1
            wrong[f"{want}->{got}"] += 1
    out["by_rule"] = dict(sorted(by_rule.items()))
    out["by_type"] = dict(sorted(by_type.items()))
    out["by_founder"] = {k: dict(sorted(v.items())) for k, v in sorted(by_founder.items())}
    out["wrong_by_kind"] = dict(wrong.most_common())
    out["unattributed"] = unattributed
    return out


def counts(d: dict) -> dict[str, int]:
    t = d["by_type"]
    return {k: t.get(k, 0) for k in ("Intestine", "Muscle", "Reproductive")} | {"decided": d["decided"]}


def founder_types(body: Body, targets: dict[str, str], names) -> dict[str, str]:  # noqa: ANN001
    """Each founder's identity: the last founder decision that fired on it, or the type it holds."""
    out = {}
    for n in names:
        c = body.cells.get(n)
        if c is None:
            out[n] = "absent"
            continue
        fired = [d for d in c.fired if d in FOUNDER_DECISIONS]
        out[n] = targets[fired[-1]] if fired else c.cell_type
    return out


# ---- registration and run ---------------------------------------------------------------------------


def sha256(path: str | Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def manifest(entries: list[dict[str, Any]], parameters: dict[str, Any]) -> dict[str, Any]:
    m = pv.manifest(
        entries,
        parameters,
        exclusions=[
            "terminal cells of AB, MS and C: no founder-tissue rule names them, so they are counted as "
            "not decided, never as wrong",
            "post-embryonic cells: the 555 are those born before 800 min",
        ],
        partitions=pv.NO_SPLIT,
        atlas=True,
    )
    m["sources"].append(CSHAPER)
    m["code_cleanliness"] = mf.code_cleanliness(__file__, OWN_CODE)
    return m


def register() -> None:
    entries = pv.inputs(atlas=False, programs=[ORG])
    payload = {
        "result": REGISTRATION_NAME,
        "date": date.today().isoformat(),
        "registration": REGISTRATION,
        "frozen": {p: sha256(p) for p in FROZEN},
        "script_sha256_at_registration": sha256(__file__),
    }
    print("saved", save_result(REGISTRATION_NAME, payload, manifest=manifest(entries, {"mode": "register"})))


def committed(path: Path) -> bool:
    tracked = subprocess.run(["git", "ls-files", "--error-unmatch", str(path)], capture_output=True)
    clean = subprocess.run(["git", "diff", "--quiet", "HEAD", "--", str(path)])
    return tracked.returncode == 0 and clean.returncode == 0


def grow(module, until: float = HATCH, **kw) -> Body:  # noqa: ANN001
    return Body(module, **kw).run(until=until)


def run() -> int:
    reg_path = RESULTS_DIR / f"{REGISTRATION_NAME}.json"
    if not reg_path.exists() or not committed(reg_path):
        raise SystemExit(f"{reg_path} is not committed: register first, commit, then run")
    reg = json.loads(reg_path.read_text())
    if reg["registration"] != json.loads(json.dumps(REGISTRATION)):
        raise SystemExit("REGISTRATION differs from the committed registration: refused")
    moved = sorted(p for p, h in reg["frozen"].items() if sha256(p) != h)
    if moved:
        raise SystemExit(f"registered files changed since registration: {moved}")
    entries = pv.inputs(atlas=False, programs=[ORG], extra=[reg_path])

    ref = ReferenceLineage.load()
    cells = terminal_cells(ref)
    module = parse_file(PROGRAM)
    targets = {d.id: d.to for d in module.decisions if d.action == "differentiate"}
    v1 = lookup_fates_present(module, cells)
    if v1:
        raise SystemExit(f"V1 void: lookup fates present: {v1[:5]}")

    out: dict[str, Any] = {
        "result": RESULT_NAME,
        "date": date.today().isoformat(),
        "registration": str(reg_path),
        "script_changed_since_registration": sha256(__file__) != reg["script_sha256_at_registration"],
    }

    # primary arm
    body = grow(module, means=True)
    diff = compare(body, ref, until=HATCH)
    primary = decided(body, cells, ref, targets)
    ref_by_800 = sum(1 for c in ref.cells.values() if c.born <= HATCH)
    void = []
    if diff.missing or diff.parent_mismatch or diff.matched != ref_by_800:
        void.append(f"V2: missing {len(diff.missing)}, parent mismatches {len(diff.parent_mismatch)}")
    if diff.fates_correct != primary["right"]:
        void.append(f"V3: compare {diff.fates_correct} against {primary['right']}")
    if primary["unattributed"] or set(primary["by_rule"]) - set(TISSUE_RULES):
        void.append(f"V4: unattributed {primary['unattributed'][:5]}, rules {primary['by_rule']}")
    out["primary"] = primary | {"compare": diff.to_dict()}
    founders = founder_types(body, targets, EXPECTED_FOUNDERS)
    out["founders_measured"] = {
        "types": founders,
        "correct": sum(founders[n] == t for n, t in EXPECTED_FOUNDERS.items()),
        "of": len(EXPECTED_FOUNDERS),
    }

    # seeds
    per_seed = {}
    for s in SEEDS:
        d = decided(grow(module, seed=s), cells, ref, targets)
        per_seed[s] = {"decided": d["decided"], "right": d["right"], "wrong": d["wrong"]}
    out["seeds"] = {
        "per_seed": per_seed,
        "decided_range": [
            min(v["decided"] for v in per_seed.values()),
            max(v["decided"] for v in per_seed.values()),
        ],
        "right_range": [
            min(v["right"] for v in per_seed.values()),
            max(v["right"] for v in per_seed.values()),
        ],
    }

    # other tables, no table, rearranged
    tables = {}
    for label, path in (("curated", CURATED), ("detected_sheet1", DETECTED)):
        m = parse_file(PROGRAM)
        m.organism.contacts = str(ROOT / path)
        b = grow(m, means=True)
        d = decided(b, cells, ref, targets)
        f = founder_types(b, targets, EXPECTED_FOUNDERS)
        tables[label] = {
            "decided": d["decided"],
            "right": d["right"],
            "wrong": d["wrong"],
            "founders": f,
            "founders_correct": sum(f[n] == t for n, t in EXPECTED_FOUNDERS.items()),
        }
    out["other_tables"] = tables
    m = parse_file(PROGRAM)
    m.organism.contacts = ""
    b = grow(m, means=True)
    d = decided(b, cells, ref, targets)
    out["no_table"] = {
        "counts": counts(d),
        "right": d["right"],
        "founders": founder_types(b, targets, EXPECTED_FOUNDERS),
    }
    swap = {"ABa": "ABp", "ABp": "ABa"}
    rows = []
    for line in MEASURED.read_text().splitlines():
        r = line.split("#", 1)[0].split()
        if len(r) >= 3:
            rows.append("\t".join([r[0], swap.get(r[1], r[1]), swap.get(r[2], r[2]), *r[3:]]))
    with tempfile.TemporaryDirectory() as tmp:  # outside data/: a scratch table, not an input
        swapped_path = Path(tmp) / "rearranged.tsv"
        swapped_path.write_text("\n".join(rows) + "\n")
        m = parse_file(PROGRAM)
        m.organism.contacts = str(swapped_path)
        f = founder_types(grow(m, until=FOUNDER_UNTIL, means=True), targets, ("ABa", "ABp"))
    out["rearranged"] = {
        "what": "ABa and ABp swapped in the measured table alone; the rules are unchanged",
        "founders": f,
        "fates_follow_position": f == {"ABa": "ABpPrecursor", "ABp": "ABaPrecursor"},
    }

    # knockouts
    kos = {}
    for factor, wanted in KNOCKOUT_FOUNDERS.items():
        b = grow(module, means=True, knockouts={factor})
        f = founder_types(b, targets, sorted(set(wanted) | {"P2"}))
        d = decided(b, cells, ref, targets)
        c = counts(d)
        test = PUBLISHED_TEST.get(factor)
        kos[factor] = {
            "founders": f,
            "founder_phenotype_reproduced": all(f.get(n) == t for n, t in wanted.items()),
            "terminal": c,
            "terminal_by_founder": d["by_founder"],
            "terminal_as_registered": c == KNOCKOUT_TERMINAL[factor],
            "published_terminal_phenotype_reproduced": test(c, d["by_founder"]) if test else None,
            "reading": REGISTRATION["published_terminal_phenotypes"][factor],
        }
    out["knockouts"] = kos

    # the registered predictions, each checked
    p = {
        "P1 primary": counts(primary) == WT_COUNTS
        and primary["right"] == 56
        and primary["wrong"] == 0
        and primary["by_rule"] == {"tissue_D": 20, "tissue_E": 34, "tissue_germline": 2},
        "P2 seeds": all(v == {"decided": 56, "right": 56, "wrong": 0} for v in per_seed.values()),
        "P3 founders": out["founders_measured"]["correct"] == 8,
        "P4 other tables": all(
            t["founders_correct"] == 8 and (t["decided"], t["right"], t["wrong"]) == (56, 56, 0)
            for t in tables.values()
        ),
        "P5 no table": out["no_table"]["counts"]
        == {"Intestine": 0, "Muscle": 20, "Reproductive": 2, "decided": 22}
        and out["no_table"]["right"] == 22
        and out["no_table"]["founders"]["ABp"] == "ABaPrecursor"
        and out["no_table"]["founders"]["E"] == "MSPrecursor",
        "P6 rearranged": out["rearranged"]["fates_follow_position"],
        "P7 knockouts, founder level": sum(k["founder_phenotype_reproduced"] for k in kos.values()) == 8,
        "P8 knockouts, terminal level": all(k["terminal_as_registered"] for k in kos.values()),
    }
    out["predictions_held"] = p
    out["void"] = void
    out["headline"] = (
        "VOID: " + "; ".join(void)
        if void
        else REGISTRATION["wording"]["P1 holds"]
        if p["P1 primary"]
        else f"{primary['decided']} of 555 decided without the lookup, {primary['right']} right, "
        f"{primary['wrong']} wrong, by rule {primary['by_rule']}"
    )
    params = {"until_min": HATCH, "seeds": SEEDS, "founder_until_min": FOUNDER_UNTIL, "timers": "means"}
    print("saved", save_result(RESULT_NAME, out, manifest=manifest(entries, params)))
    print(json.dumps({k: out[k] for k in ("headline", "predictions_held", "void")}, indent=1))
    return 1 if void else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--write-tree", action="store_true", help="derive lineage_embryo_tree.bio")
    ap.add_argument("--register", action="store_true", help="write the registration")
    a = ap.parse_args()
    if a.write_tree:
        TREE.write_text(tree_text(LINEAGE.read_text()))
        print("wrote", TREE)
        return 0
    if a.register:
        register()
        return 0
    return run()


if __name__ == "__main__":
    sys.exit(main())
