"""Digital Development (Du et al. 2014): the observed founder transformations per knockout, scored against
the program's own knockouts."""

from pathlib import Path

import pytest

from genomeos.lang import parse_file
from genomeos.organism.digital_development import GENE_TO_FACTOR, IDENTITY_TYPE, compare
from genomeos.results import load_result


def test_compare_scores_reproduced_missed_and_not_modelled():
    m = parse_file("data/organisms/celegans/embryo.bio")
    table = {
        "pop-1": {"MS": "E", "ABpla": "ABpra"},  # MS adopts E (modelled); an 8-cell change (not modelled)
        "apx-1": {"ABp": "ABa"},
        "pie-1": {"P2": "EMS", "ABp": "ABa"},
        "skn-1": {"EMS": "C"},  # observed at the parent: the program types the daughters
        "pal-1": {"C": "D"},  # deliberately wrong: the program says C loses its fate, not that it becomes D
        "zyg-11": {"AB": "P1"},  # not a factor the program has: skipped
    }
    r = compare(m, table, until=300)
    by = {g["gene"]: g for g in r["genes"]}
    assert set(by) == {"pop-1", "apx-1", "pie-1", "skn-1", "pal-1"} and r["genes_in_table"] == 6
    assert [o["status"] for o in by["pop-1"]["observed"]] == ["reproduced", "not modelled"]
    assert [o["status"] for o in by["apx-1"]["observed"]] == ["reproduced"]
    assert [o["status"] for o in by["pie-1"]["observed"]] == ["reproduced", "reproduced"]
    assert [o["status"] for o in by["skn-1"]["observed"]] == ["reproduced"]
    assert [o["status"] for o in by["pal-1"]["observed"]] == ["missed"]
    assert (
        r["reproduced"] == 5
        and r["missed"] == 1
        and r["not_modelled"] == 1
        and r["unmodelled_genes"] == ["zyg-11"]
    )
    assert {p["mutant"] for p in by["pop-1"]["predicted"]} == {"EPrecursor"}
    assert (
        set(IDENTITY_TYPE) == {"ABa", "ABp", "EMS", "MS", "E", "C", "D"}
        and GENE_TO_FACTOR["mom-2"] == "MOM-2"
    )


@pytest.mark.skipif(
    not Path("data/results/celegans_digital_development.json").exists(),
    reason="run `genomeos data distil --only celegans_digital_development`",
)
def test_distilled_digital_development_check_holds():
    r = load_result("celegans_digital_development")
    assert r["genes_in_table"] >= 70 and r["genes_modelled"] >= 7
    by = {g["gene"]: g for g in r["genes"]}
    for gene, cell, identity in (
        ("pop-1", "MS", "E"),
        ("apx-1", "ABp", "ABa"),
        ("pie-1", "P2", "EMS"),
        ("mom-2", "E", "MS"),
    ):
        status = next(
            o["status"] for o in by[gene]["observed"] if o["cell"] == cell and o["adopts"] == identity
        )
        assert status == "reproduced", (gene, cell, identity)
    # the misses are par-2 and par-3, whose polarity rules are not written
    assert r["rate_over_modelled"] >= 0.6
