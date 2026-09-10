import os

import pytest

from genomeos.molecules import parse_pdb

PDB = """ATOM      1  N   MET A   1      -8.901   4.127  -0.555  1.00 95.20           N
ATOM      2  CA  MET A   1      -8.608   3.135  -1.618  1.00 95.20           C
ATOM      3  CA  ALA A   2      -5.123   1.024  -0.100  1.00 40.00           C
ATOM      4  CA  GLY A   3      -2.000   0.000   0.000  1.00 80.00           C
"""


def test_parse_pdb_reads_ca_and_plddt():
    s = parse_pdb(PDB, "AF-TEST-F1-model_v6", "TEST")
    assert s.length == 3 and "".join(s.residues) == "MAG"
    assert s.plddt == [95.2, 40.0, 80.0] and abs(s.mean_plddt() - 71.73) < 0.01
    assert abs(s.confident_fraction() - 2 / 3) < 1e-9 and s.radius_of_gyration() > 0
    d = s.to_dict()
    assert d["confidence"] == 0.72 and "predicted" in d["evidence"]


@pytest.mark.skipif(
    not os.environ.get("GENOMEOS_NET"), reason="set GENOMEOS_NET=1 for live UniProt/AlphaFold calls"
)
def test_live_app_report():
    from genomeos.molecules import protein_report

    r = protein_report("APP", structure=True)
    assert r["uniprot"]["accession"] == "P05067" and r["uniprot"]["length"] == 770
    assert r["structure"]["length"] == 770 and r["structure"]["sequence_matches_uniprot"]
