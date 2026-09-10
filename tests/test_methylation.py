import gzip

from genomeos.twin.clocks import Clock
from genomeos.twin.methylation import load_probes, summarise


def test_probe_map_is_complete():
    probes = load_probes()
    assert len(probes) == 418
    assert set(Clock.horvath().coefficients) <= set(probes)
    assert all(chrom.startswith("chr") and beg > 0 for chrom, beg in probes.values())


def test_summarise_pools_strands_and_files(tmp_path):
    probes = load_probes()
    pid, (chrom, beg) = next(iter(probes.items()))

    def row(pos, cov, pct, code="m"):
        tail = "\t".join(["0"] * 7)
        head = f"{chrom}\t{pos}\t{pos + 1}\t{code}\t{cov}\t.\t{pos}\t{pos + 1}\t255,0,0"
        return f"{head}\t{cov}\t{pct:.2f}\t{tail}\n"

    f1 = tmp_path / "hap1.bedmethyl.gz"
    f2 = tmp_path / "hap2.bedmethyl.gz"
    with gzip.open(f1, "wt") as fh:
        fh.write(row(beg, 10, 80.0) + row(beg + 1, 10, 60.0))  # both strands
        fh.write(
            f"{chrom}\t{beg}\t{beg + 1}\th\t3\t.\t{beg}\t{beg + 1}\t255,0,0\t3\t10.00\t0\t0\t0\t0\t0\t0\t0\n"
        )  # 5hmC ignored
    with gzip.open(f2, "wt") as fh:
        fh.write(row(beg, 20, 100.0))
    s = summarise([str(f1), str(f2)], probes, min_coverage=5)
    assert s.probes_covered == 1 and s.rows_scanned == 3
    assert abs(s.betas[pid] - (10 * 0.8 + 10 * 0.6 + 20 * 1.0) / 40) < 1e-9
    r = Clock.horvath().predict(s.betas)
    assert r.missing == 352 and 0 < r.age < 200
