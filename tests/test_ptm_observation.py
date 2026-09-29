# SPDX-License-Identifier: AGPL-3.0-or-later
"""The phosphosite observation join: accession, position and residue, phospho class only."""

from genomeos.molecules import ptm


def _defn(acc, seq, feats):
    return {
        "id": f"UniProt:{acc}",
        "gene": acc,
        "sections": {
            "identity": {"items": {"accession": acc, "sequence": seq}},
            "modifications": {"items": feats},
        },
    }


def _mod(desc, pos, kind="Modified residue"):
    return {"type": kind, "description": desc, "start": pos, "end": pos}


def test_join_counts_only_phospho_with_matching_residue():
    seq = "MSTYCKS"
    d = _defn(
        "P1",
        seq,
        [
            _mod("Phosphoserine", 2),  # S2: in reference, residue agrees -> joined
            _mod("Phosphothreonine", 3),  # T3: reference says S -> disagreement, not joined
            _mod("Phosphotyrosine", 4),  # Y4: not in reference
            _mod("Disulfide", 5, "Disulfide bond"),  # C5: not phospho, never joined
            _mod("N6-acetyllysine", 6),  # K6: other class on a reference position
        ],
    )
    ref = {("P1", 2): ("S", 7, 40), ("P1", 3): ("S", 2, 5), ("P1", 5): ("C", 1, 1), ("P1", 6): ("K", 1, 1)}
    j = ptm.join_observations(ref, [d])
    assert j["curated_sites"] == 5
    assert j["curated_phospho_sites"] == 3
    assert j["accession_position_matches"] == 2
    assert j["residue_disagreements"] == 1
    assert j["joined"] == 1
    assert j["sites"] == {"P1": [[2, "S", 7, 40]]}
    assert j["other_class_on_reference_position"] == {"disulfide": 1, "acetyl": 1}
    assert j["off_acceptor"] == []


def test_observed_in_names_the_field_and_its_limit():
    obs = {"sites": {"P1": [[2, "S", 7, 40]]}}
    rec = ptm.observed_in("P1", 2, obs)
    assert rec["observed_in_cell_types_or_tissues"] == 7
    assert rec["spectral_count"] == 40
    assert "not occupancy" in rec["meaning"]
    assert ptm.observed_in("P1", 3, obs) is None
    assert ptm.observed_in("P1", 2, None) is None
    assert not any(k in rec for k in ("occupancy", "modified", "active"))
