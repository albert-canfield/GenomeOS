# SPDX-License-Identifier: AGPL-3.0-or-later
"""Two people's protein-changing variants as a code review: normalised keys, evidence, genes ranked."""

import json

from genomeos.genome import diff


def _person(root, name, variants, scores=None, clinvar=None):
    d = root / name
    d.mkdir(parents=True)
    (d / "coding_variants.json").write_text(json.dumps(variants))
    if scores:
        (d / "alphamissense.json").write_text(json.dumps({"scores": scores, "asked": list(scores)}))
    if clinvar:
        (d / "clinvar_screen.json").write_text(json.dumps({"hits": clinvar}))


def _v(gene, chrom, pos, ref, alt, cons, hgvs, gt):
    hom = gt.replace("|", "/").split("/").count("1") >= 2
    return {
        "gene": gene,
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "consequence": cons,
        "hgvs_p": hgvs,
        "genotype": gt,
        "zygosity": "homozygous" if hom else "heterozygous",
        "transcript": "T1",
    }


def test_diff_normalises_attaches_evidence_and_ranks(tmp_path):
    a = [
        _v("APP", "chr21", 100, "TGG", "TGGG", "frameshift", "p.Gly5fs", "0/1"),
        _v("SOD1", "chr21", 200, "C", "T", "missense", "p.Arg2Trp", "1/1"),
        _v("DYRK1A", "chr21", 300, "G", "A", "missense", "p.Ala9Thr", "0/1"),
    ]
    b = [
        _v(
            "APP", "chr21", 100, "T", "TG", "frameshift", "p.Gly5fs", "0/1"
        ),  # the same insertion, written differently
        _v("DYRK1A", "chr21", 300, "G", "A", "missense", "p.Ala9Thr", "0/1"),
        _v("KCNE1", "chr21", 400, "A", "G", "nonsense", "p.Trp3Ter", "0/1"),
    ]
    _person(
        tmp_path,
        "A",
        a,
        scores={
            "chr21:200:C:T": [
                {"transcript": "T1", "score": 0.91, "class": "likely_pathogenic", "protein_variant": "R2W"}
            ]
        },
        clinvar=[
            {
                "chrom": "chr21",
                "pos": 200,
                "ref": "C",
                "alt": "T",
                "significance": "Pathogenic",
                "conditions": "ALS",
            }
        ],
    )
    _person(tmp_path, "B", b)
    d = diff.genome_diff("A", "B", tmp_path)
    assert d["variants_a"] == 3 and d["variants_b"] == 3 and d["shared"] == 2
    assert d["only_a"] == 1 and d["only_b"] == 1 and d["genes_differing"] == 2
    genes = {g["gene"]: g for g in d["genes"]}
    assert "APP" not in genes  # one insertion written two ways is one variant, so APP does not differ
    sod1 = genes["SOD1"]["only_a"][0]
    assert (
        sod1["predicted"]["class"] == "likely_pathogenic" and sod1["clinvar"]["significance"] == "Pathogenic"
    )
    assert d["likely_pathogenic_only_a"] == 1 and d["clinvar_only_a"] == 1 and d["truncating_only_b"] == 1
    assert d["genes"][0]["gene"] == "SOD1"  # ClinVar first, then truncations
    md = diff.render(d)
    assert "+ p.Arg2Trp (missense, homozygous, 1/1)  [likely pathogenic 0.91; ClinVar Pathogenic: ALS]" in md
    assert "- p.Trp3Ter (nonsense, heterozygous, 0/1)" in md and "A against B" in md


def test_diff_against_the_reference_is_the_persons_own_edits(tmp_path):
    _person(tmp_path, "A", [_v("APP", "chr21", 100, "C", "T", "missense", "p.Arg1Trp", "0/1")])
    d = diff.genome_diff("A", "reference", tmp_path)
    assert d["variants_b"] == 0 and d["only_a"] == 1 and d["genes_differing"] == 1
    assert "against the reference" in diff.render(d)
