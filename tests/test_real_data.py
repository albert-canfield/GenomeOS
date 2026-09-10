"""Tests against real human sequence. Skipped unless the files have been
fetched with scripts/fetch_reference.sh (they are gitignored)."""

from pathlib import Path

import pytest

from genomeos.genome import read_fasta
from genomeos.runtime import STANDARD_START_CODONS, VERTEBRATE_MITOCHONDRIAL_CODE, find_orfs

CHRM = Path("data/reference/chrM.fa.gz")


@pytest.mark.skipif(not CHRM.exists(), reason="run scripts/fetch_reference.sh chrM")
def test_human_mitochondrial_genome_translates_to_known_proteins():
    seq = read_fasta(CHRM)["chrM"]
    assert len(seq) == 16569  # rCRS length
    # AUG-only starts reproduce the exact published lengths of the AUG-initiated genes;
    # the mitochondrial initiators (AUU/AUA) are covered in tests/test_annotation.py
    orfs = {
        o.locus.start: o
        for o in find_orfs(
            seq, "chrM", min_aa=200, table=VERTEBRATE_MITOCHONDRIAL_CODE, start_codons=STANDARD_START_CODONS
        )
    }
    # MT-CO1 (cytochrome c oxidase I): 513 aa, ATG at 5903 (0-based)
    assert orfs[5903].length_aa == 513
    assert orfs[5903].protein.startswith("MFADRWLFSTNHKDIGTLYLLFGAWAGVLGTALSLLIRAELGQPGNLLGNDHIYNVIVTA")
    # MT-CO2: 227 aa at 7585; MT-ATP6: 226 aa at 8526
    assert orfs[7585].length_aa == 227
    assert orfs[8526].length_aa == 226
