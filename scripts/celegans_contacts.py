# SPDX-License-Identifier: AGPL-3.0-or-later
"""The worm's founders decided by contact instead of by name (area E, BioLang v0.4 §7.4).

embryo_contacts.bio is embryo_factors.bio with founders_contacts.bio in place of founders.bio and a
contact table on the organism block. Its three founder signals name no sender and no receiver: each
receiver reads the ligand summed over the cells touching it, and who touches whom comes from
data/organisms/celegans/contacts_embryo.tsv. This script measures what that buys and what it costs:

1. the seven founder fates, against Sulston;
2. the rearrangement Priess & Thomson 1987 and Hutter & Schnabel 1994 did — ABa put where ABp is —
   which the named-sender program cannot express at all: the fates must follow the contact;
3. the eight published knockouts, which must still come out as they did with named senders;
4. the whole embryo against the reference lineage, cell by cell, before and after: the number area E
   is judged on must not move because the mechanism underneath it changed.

    uv run python scripts/celegans_contacts.py
"""

import tempfile
import time
from pathlib import Path

from genomeos.lang import parse_file
from genomeos.organism import REFERENCES, ReferenceLineage
from genomeos.organism.diff import compare
from genomeos.results import save_result
from genomeos.runtime.body import Body

NAMED = "data/organisms/celegans/embryo_factors.bio"
CONTACTS = "data/organisms/celegans/embryo_contacts.bio"
TABLE = "data/organisms/celegans/contacts_embryo.tsv"
FOUNDERS = ["ABa", "ABp", "EMS", "P2", "MS", "E", "C", "D", "P4"]
UNTIL = 800.0

# what Sulston 1983 and the genetics say each founder is, as the program's own type names
EXPECTED = {
    "ABa": "ABaPrecursor",
    "ABp": "ABpPrecursor",
    "EMS": "EMSPolarised",
    "MS": "MSPrecursor",
    "E": "EPrecursor",
    "C": "CPrecursor",
    "D": "DPrecursor",
    "P4": "GermCell",
}

# the published phenotype of each knockout, as the founder types it should produce
KNOCKOUTS = {
    "APX-1": ("ABp becomes ABa-like: no Delta reaches GLP-1", {"ABa": "ABaPrecursor", "ABp": "ABaPrecursor"}),
    "GLP-1": ("ABp becomes ABa-like: the receptor is gone", {"ABa": "ABaPrecursor", "ABp": "ABaPrecursor"}),
    "MOM-2": ("E takes the MS fate: no Wnt at the EMS/P2 contact", {"MS": "MSPrecursor", "E": "MSPrecursor"}),
    "MOM-5": ("E takes the MS fate: the receptor is gone", {"MS": "MSPrecursor", "E": "MSPrecursor"}),
    "PIE-1": (
        "P2 takes the EMS fate, and with the P2 identity lost ABp becomes ABa-like",
        {"P2": "EMSPrecursor", "ABp": "ABaPrecursor"},
    ),
    "SKN-1": ("EMS daughters take the C fate through PAL-1", {"MS": "CPrecursor", "E": "CPrecursor"}),
    "POP-1": ("MS takes the E fate", {"MS": "EPrecursor", "E": "EPrecursor"}),
    "PAL-1": ("C and D lose their fates", {"C": "Blastomere", "D": "Blastomere"}),
}


def grow(module, until=UNTIL, **kw):
    return Body(module, **kw).run(until=until)


def founders(body):
    return {n: body.cells[n].cell_type for n in FOUNDERS if n in body.cells}


def rearranged_table(path, out):
    """ABa where ABp is and ABp where ABa is: the blastomere rearrangement, as a contact table."""
    swap = {"ABa": "ABp", "ABp": "ABa"}
    lines = []
    for line in Path(path).read_text().splitlines():
        body, _, comment = line.partition("#")
        row = body.split()
        if len(row) >= 3:
            row[1], row[2] = swap.get(row[1], row[1]), swap.get(row[2], row[2])
            body = "\t".join(row)
        lines.append(body + ("#" + comment if comment else ""))
    Path(out).write_text("\n".join(lines) + "\n")
    return out


REFERENCE = ReferenceLineage.load(REFERENCES["celegans"])


def score(module, **kw):
    """Cells, fates, deaths and timing against the reference lineage, as `genomeos grow --compare`."""
    body = grow(module, **kw)
    return body, compare(body, REFERENCE, until=UNTIL).to_dict()


t0 = time.time()
named = parse_file(NAMED)
contacts = parse_file(CONTACTS)

# 1. the founders, and the same run scored against the reference lineage
named_body, named_diff = score(named)
body, diff = score(contacts)
result = {
    "programs": {"named_senders": NAMED, "contacts": CONTACTS, "table": TABLE},
    "founders": {
        "with_contacts": founders(body),
        "with_named_senders": founders(named_body),
        "expected": EXPECTED,
        "correct": sum(1 for n, t in EXPECTED.items() if founders(body).get(n) == t),
        "of": len(EXPECTED),
    },
}

# 2. the rearrangement: nothing about the rules changes, only which cells touch
with tempfile.TemporaryDirectory() as tmp:
    swapped = parse_file(CONTACTS)
    swapped.organism.contacts = rearranged_table(TABLE, Path(tmp) / "rearranged.tsv")
    sw = grow(swapped, until=250)
    # and the control: no table at all, so no cell touches any other
    blind = parse_file(CONTACTS)
    blind.organism.contacts = ""
    bl = grow(blind, until=250)
result["rearranged"] = {
    "what": "ABa put where ABp is and ABp where ABa is, in the contact table alone",
    "founders": founders(sw),
    "fates_follow_position": sw.cells["ABa"].cell_type == "ABpPrecursor"
    and sw.cells["ABp"].cell_type == "ABaPrecursor",
    "published": "Priess & Thomson 1987, Cell 48:241; Hutter & Schnabel 1994, Development 120:2051",
}
result["no_contact_table"] = {
    "what": "the same program with no table, so nothing touches anything: a control on the table",
    "founders": founders(bl),
    "nothing_induced": bl.cells["ABp"].cell_type == "ABaPrecursor"
    and bl.cells["E"].cell_type == "MSPrecursor",
}

# 3. the eight published knockouts
kos = {}
for factor, (expect, wanted) in KNOCKOUTS.items():
    got = founders(grow(contacts, until=250, knockouts={factor}))
    kos[factor] = {
        "expect": expect,
        "founders": got,
        "ok": all(got.get(n) == t for n, t in wanted.items()),
    }
result["knockouts"] = kos
result["knockouts_reproduced"] = sum(1 for v in kos.values() if v["ok"])

# 4. the whole embryo, cell by cell, before and after
fa = {n: (c.cell_type, c.terminal_name, round(c.born, 6)) for n, c in named_body.cells.items()}
fb = {n: (c.cell_type, c.terminal_name, round(c.born, 6)) for n, c in body.cells.items()}
differing = sorted(n for n in set(fa) & set(fb) if fa[n] != fb[n])
result["against_sulston"] = {
    "named_senders": named_diff,
    "contacts": diff,
    "cells_only_in_one": sorted(set(fa) ^ set(fb)),
    "cells_that_differ": differing,
    "difference": {n: {"named": fa[n], "contacts": fb[n]} for n in differing},
}
result["seconds"] = round(time.time() - t0, 1)
result["pass"] = bool(
    result["founders"]["correct"] == result["founders"]["of"]
    and result["rearranged"]["fates_follow_position"]
    and result["no_contact_table"]["nothing_induced"]
    and result["knockouts_reproduced"] == len(KNOCKOUTS)
    and diff["fates_correct"] >= named_diff["fates_correct"]
    and diff["matched"] == named_diff["matched"]
    and not result["against_sulston"]["cells_only_in_one"]
)
print(save_result("celegans_contacts", result))
print("founders", result["founders"]["correct"], "/", result["founders"]["of"])
print("rearranged", result["rearranged"]["founders"])
print("knockouts", result["knockouts_reproduced"], "/", len(KNOCKOUTS))
print("fates", named_diff["fates_correct"], "->", diff["fates_correct"], "of", diff["fates_checked"])
print("cells that differ:", differing)
print("pass", result["pass"])
