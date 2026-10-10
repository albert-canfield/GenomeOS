# SPDX-License-Identifier: AGPL-3.0-or-later
"""Area E's fates decided without the lookup: the measured contact table, the lineage with its fates
taken out, and the instrument that counts (scripts/celegans_nolookup.py).

The instrument is tested on a toy before it is trusted on the worm, because the number it reports is a
count of fates the lookup did NOT write, and an instrument that let one lookup fate through, or
credited a type a cell merely inherited without saying from where, would report a larger number and
be wrong.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

from genomeos.lang import parse, parse_file
from genomeos.organism.reference import ReferenceLineage
from genomeos.runtime.body import Body, read_contacts
from scripts.celegans_nolookup import (
    LINEAGE,
    PROGRAM,
    TREE,
    attribute,
    counts,
    decided,
    lookup_fates_present,
    terminal_cells,
    tree_text,
)
from tests.committed_data import committed

ORG = Path("data/organisms/celegans")
MEASURED = ORG / "contacts_embryo_cshaper.tsv"
DETECTED = ORG / "contacts_embryo_cshaper_detected.tsv"


def _partners(path: Path, cell: str) -> set[str]:
    (snap,) = read_contacts(str(path))  # one snapshot: the source states no times
    assert snap[0] == 0.0
    return set(snap[1].get(cell, {}))


# ---- the measured table ------------------------------------------------------------------------------


def test_the_measured_tables_have_the_published_pair_counts_and_no_unknown_name():
    ref = ReferenceLineage.load()
    for path, pairs in ((MEASURED, 479), (DETECTED, 3011)):
        rows = [r.split("\t") for r in path.read_text().splitlines() if r and not r.startswith("#")]
        assert len(rows) == pairs
        assert all(r[0] == "0" and r[3] == "1.0" for r in rows)  # topology only, no time, no area
        assert {c for r in rows for c in r[1:3]} <= set(ref.cells)
        assert "sha256 519c34349af8e8f6cbd675c665f9a880c99fe1ad482b469067f4a6078a855c99" in path.read_text()


def test_aba_never_touches_p2_in_any_measured_reading():
    """The one geometric fact that separates ABp from ABa (Priess & Thomson 1987), now measured: not one
    of 17 embryos shows ABa touching P2 at any time, at any contact area (Sheet1)."""
    for path in (MEASURED, DETECTED):
        assert _partners(path, "ABa") == {"ABp", "EMS"}
        assert {"ABp", "EMS"} <= _partners(path, "P2")
    assert "P2" not in _partners(MEASURED, "MS")  # MS is born away from the Wnt source
    assert {"ABalp", "ABara"} <= _partners(MEASURED, "MS")  # the second Notch interaction's receivers


# ---- the lineage without its fates -------------------------------------------------------------------


def test_the_tree_is_exactly_the_lineage_without_its_fate_decisions():
    lineage = LINEAGE.read_text()
    tree = TREE.read_text()
    assert tree == tree_text(lineage)
    assert not any(line.startswith("decision fate_") for line in tree.splitlines())
    kinds = ("decision div_", "decision die_")
    for k in kinds:
        assert sum(x.startswith(k) for x in tree.splitlines()) == sum(
            x.startswith(k) for x in lineage.splitlines()
        )
    assert sum(x.startswith("decision fate_") for x in lineage.splitlines()) == 555


def test_the_program_holds_no_lookup_fate_for_any_of_the_555():
    cells = terminal_cells(ReferenceLineage.load())
    assert len(cells) == 555
    assert lookup_fates_present(parse_file(PROGRAM), cells) == []
    # and the instrument does see one when it is there
    assert lookup_fates_present(parse_file(ORG / "embryo_contacts.bio"), cells)


# ---- the instrument, on a toy ------------------------------------------------------------------------

TOY = """
module toy.nolookup
import bio.std.development
organism O { root: Z; cell_type: Zygote }
stage S { from: 0 min }
cell_type EPrecursor { parent: Blastomere }
cell_type Intestine { parent: Blastomere }
decision z { action: divide; when: cell = Z; daughters: A, B; after: 1 min }
decision a { action: divide; when: cell = A; daughters: Aa, Ap; after: 1 min }
decision founder_A { action: differentiate; when: cell = A; to: EPrecursor }
decision tissue_E { action: differentiate; when: cell_type = EPrecursor; to: Intestine }
"""


def _toy() -> tuple[Body, dict[str, str]]:
    m = parse(TOY)
    body = Body(m).run(until=5)
    return body, {d.id: d.to for d in m.decisions if d.action == "differentiate"}


def test_a_type_set_on_the_founder_is_attributed_to_the_rule_that_set_it():
    body, targets = _toy()
    assert body.cells["Aa"].cell_type == "Intestine" and "tissue_E" not in body.cells["Aa"].fired
    assert attribute(body, "Aa", targets) == "tissue_E"  # found on A, up the parents
    assert attribute(body, "B", targets) is None  # B was never typed


def test_the_count_separates_right_wrong_and_undecided():
    body, targets = _toy()
    ref = SimpleNamespace(cells={})  # no founders known: every cell's founder is "none"
    d = decided(body, {"Aa": "Intestine", "Ap": "Muscle", "B": "Neuron"}, ref, targets)
    assert (d["cells"], d["decided"], d["right"], d["wrong"]) == (3, 2, 1, 1)
    assert d["by_rule"] == {"tissue_E": 2} and d["unattributed"] == []
    assert d["wrong_by_kind"] == {"Muscle->Intestine": 1}
    assert d["by_founder"] == {"none": {"Intestine": 2}}


# ---- the registered result, pinned (2026-10-10, registration 6bdecd4, result 09be446) ----------------

RESULT_RELATIVE = "data/results/celegans_nolookup.json"
#: Intestine / Muscle / Reproductive / decided among the 555, per knockout, as registered and measured
KNOCKOUTS = {
    "POP-1": (111, 20, 2, 133),
    "SKN-1": (0, 20, 2, 22),
    "PIE-1": (0, 20, 0, 20),
    "APX-1": (34, 20, 2, 56),
    "GLP-1": (34, 20, 2, 56),
    "PAL-1": (34, 0, 2, 36),
    "MOM-2": (0, 20, 2, 22),
    "MOM-5": (0, 20, 2, 22),
}


def _row(c: dict) -> tuple[int, int, int, int]:
    return (c["Intestine"], c["Muscle"], c["Reproductive"], c["decided"])


@committed(RESULT_RELATIVE)
def test_the_committed_result_is_the_registered_one_and_nothing_is_void():
    d = json.loads(Path(RESULT_RELATIVE).read_text())
    p = d["primary"]
    assert (p["cells"], p["decided"], p["right"], p["wrong"]) == (555, 56, 56, 0)
    assert p["by_rule"] == {"tissue_D": 20, "tissue_E": 34, "tissue_germline": 2}
    assert p["compare"]["fates_correct"] == 56  # the Body's own instrument agrees (V3)
    assert all(d["predictions_held"].values()) and len(d["predictions_held"]) == 8
    assert d["void"] == [] and d["script_changed_since_registration"] is False
    assert {k: _row(v["terminal"]) for k, v in d["knockouts"].items()} == KNOCKOUTS
    published = {k: v["published_terminal_phenotype_reproduced"] for k, v in d["knockouts"].items()}
    assert published == {
        "APX-1": None,
        "GLP-1": None,
        "MOM-2": True,
        "MOM-5": True,
        "PIE-1": False,  # registered as not reproduced: the program's pie-1 has no intestine at all
        "SKN-1": True,
        "POP-1": True,
        "PAL-1": True,
    }
    assert d["result_manifest"]["complete"] is True


def test_the_program_still_decides_the_registered_fates_and_no_more():
    """A live re-run, so a change to the program, the founders or the runtime that moves the count fails
    here and not only at the next registered run."""
    ref = ReferenceLineage.load()
    cells = terminal_cells(ref)
    module = parse_file(PROGRAM)
    targets = {d.id: d.to for d in module.decisions if d.action == "differentiate"}
    d = decided(Body(module, means=True).run(until=800), cells, ref, targets)
    assert (d["decided"], d["right"], d["wrong"]) == (56, 56, 0)
    assert d["by_founder"] == {"D": {"Muscle": 20}, "E": {"Intestine": 34}, "P4": {"Reproductive": 2}}
    for factor, row in KNOCKOUTS.items():
        k = decided(Body(module, means=True, knockouts={factor}).run(until=800), cells, ref, targets)
        assert _row(counts(k)) == row, factor
