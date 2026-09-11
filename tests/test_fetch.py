import gzip

from genomeos.genome import annotation as ann_mod
from genomeos.genome.fetch import gencode_chrom_path


def test_default_gencode_finds_a_single_chromosome_file(tmp_path, monkeypatch):
    monkeypatch.setattr(ann_mod, "GENCODE_FULL", tmp_path / "none.gff3.gz")
    monkeypatch.setattr(ann_mod, "GENCODE_SUBSET", tmp_path / "none2.gff3.gz")
    monkeypatch.chdir(tmp_path)
    (tmp_path / "data" / "results").mkdir(parents=True)
    assert ann_mod.default_gencode({"chr22"}) is None
    p = gencode_chrom_path("chr22")
    with gzip.open(p, "wt") as fh:
        fh.write("##gff-version 3\nchr22\tHAVANA\tgene\t1\t2\t.\t+\t.\tID=g1;gene_type=lncRNA;gene_name=X\n")
    assert ann_mod.default_gencode({"chr22"}) == p
    assert ann_mod.default_gencode({"chr22", "chr21"}) is None  # one file covers one chromosome
    assert ann_mod.default_gencode({"chr3"}) is None
