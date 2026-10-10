# SPDX-License-Identifier: AGPL-3.0-or-later
"""The registration of the MaBoSS score (genomeos/molecules/maboss_reference.py), pinned as it was fixed
before any scored run: the published model files byte for byte, the printed probabilities, the
phenotype rule, the mutant rule, the tolerance, the run settings and what each engine is expected to
do. ROADMAP §5 item 13, area C."""

import hashlib
import json
from pathlib import Path

from genomeos.molecules import maboss_reference as ref

FROZEN_TOLERANCE = {"printed_phenotype_abs": 0.02, "unprinted_total_max": 0.02}
FROZEN_AGREEMENT = 0.02
FROZEN_REGISTRATION_SHA256 = "5c8b08d173cb0a5bb342a21f5ca2a4e734e5d3edccd6d719ab6cf2b094fd4cd6"


def _registration_digest() -> str:
    payload = {
        "figure": ref.FIGURE,
        "conditions": ref.CONDITIONS,
        "precedence": ref.PRECEDENCE,
        "outputs": ref.OUTPUTS,
        "expected": ref.EXPECTED,
        "settings": ref.SETTINGS,
        "engines": ref.ENGINES,
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def test_the_references_and_the_tolerance_are_the_ones_fixed_before_the_run():
    assert ref.TOLERANCE == FROZEN_TOLERANCE and ref.AGREEMENT == FROZEN_AGREEMENT
    assert _registration_digest() == FROZEN_REGISTRATION_SHA256
    assert all(abs(sum(panel.values()) - 1.0) < 1e-9 for panel in ref.FIGURE.values())
    assert sum(len(p) + 1 for p in ref.FIGURE.values()) == 16
    for kind in ("bnd", "cfg"):
        digest = hashlib.sha256(Path(ref.ROOT_FILES[kind]).read_bytes()).hexdigest()
        assert digest == ref.SOURCE[f"{kind}_sha256"], f"{kind} is not the published file"


def test_the_published_files_name_the_nodes_the_registration_reads():
    """The six outputs, the two mutated nodes and the four random inputs, read from the files as text."""
    bnd = Path(ref.ROOT_FILES["bnd"]).read_text()
    cfg = Path(ref.ROOT_FILES["cfg"]).read_text()
    nodes = {line.split()[1] for line in bnd.splitlines() if line.startswith("Node ")}
    assert len(nodes) == 32 and set(ref.OUTPUTS) <= nodes and {"p53", "NICD"} <= nodes
    commented = {
        line[2:].split(".")[0] for line in cfg.splitlines() if line.startswith("//") and ".istate" in line
    }
    assert commented == {"ECMicroenv", "DNAdamage", "GF", "TGFbeta"}
    assert all(f"{o}.is_internal=0;" in cfg for o in ref.OUTPUTS)
    assert "seed_pseudorandom = 100;" in cfg and "sample_count = 50000;" in cfg and "max_time = 50;" in cfg


def test_the_phenotype_rule_names_every_table2_state_by_its_column():
    table2 = {
        "HS": (),
        "Apoptosis": ("Apoptosis", "CellCycleArrest"),
        "EMT": ("EMT", "CellCycleArrest"),
        "Metastasis": ("Metastasis", "Migration", "Invasion", "EMT", "CellCycleArrest"),
    }
    for name, state in table2.items():
        assert ref.label(state) == name
    assert ref.label(("Invasion", "EMT", "CellCycleArrest")) == "Invasion"
    assert ref.label(("CellCycleArrest",)) == "CCA" and ref.label(("p53",)) == "HS"


def test_each_check_applies_the_bound_both_ways():
    exact = {"HS": 0.25, "Metastasis": 0.25, "EMT": 0.01, "Apoptosis": 0.49}
    assert all(c["pass"] for c in ref.score("A. Wild type", exact))
    off = {**exact, "HS": 0.25 - 0.021, "CCA": 0.021}
    checks = {c["check"]: c["pass"] for c in ref.score("A. Wild type", off)}
    assert checks["HS probability"] is False and checks["phenotypes the panel does not print, total"] is False
    inside = {"Metastasis": 0.981, "EMT": 0.019}
    assert all(c["pass"] for c in ref.score("D. p53 LoF NICD GoF", inside))
    assert ref.phenotypes({(): 0.2, ("Apoptosis", "CellCycleArrest"): 0.3, ("CellCycleArrest",): 0.5}) == {
        "HS": 0.2,
        "Apoptosis": 0.3,
        "CCA": 0.5,
    }
    assert abs(ref.total_variation({("A",): 0.5, (): 0.5}, {("A",): 0.4, (): 0.6}) - 0.1) < 1e-12
