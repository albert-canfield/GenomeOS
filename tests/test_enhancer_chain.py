"""The chromosome chain decides completion from the scoring script's own result flag."""

import importlib.util
import json
from pathlib import Path


def _chain():
    spec = importlib.util.spec_from_file_location("chain", Path("scripts/enhancer_targets_all_chain.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_complete_and_scored_read_the_result_flag(tmp_path, monkeypatch):
    chain = _chain()
    monkeypatch.setattr(chain, "RESULTS", tmp_path)
    assert chain.complete("chr1") is False and chain.scored("chr1") is None
    (tmp_path / "enhancer_targets_all_chr1.json").write_text(
        json.dumps({"scored": 200, "elements_total": 900, "complete": False})
    )
    assert chain.complete("chr1") is False and chain.scored("chr1") == (200, 900)
    (tmp_path / "enhancer_targets_all_chr1.json").write_text(
        json.dumps({"scored": 900, "elements_total": 900, "complete": True})
    )
    assert chain.complete("chr1") is True
    (tmp_path / "enhancer_targets_all_chr2.json").write_text("{not json")
    assert chain.complete("chr2") is False and chain.scored("chr2") is None
    # smallest first, every chromosome but chr21 (already running) and chrM, each once
    assert chain.ORDER[:3] == ["chr22", "chrY", "chr19"] and chain.ORDER[-1] == "chr1"
    expected = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
    assert set(chain.ORDER) == expected - {"chr21"} and len(chain.ORDER) == len(set(chain.ORDER))


def test_committable_refuses_the_per_element_table(tmp_path, monkeypatch):
    chain = _chain()
    monkeypatch.setattr(chain, "RESULTS", tmp_path)
    p = tmp_path / "enhancer_targets_all_chr22.json"
    p.write_text(json.dumps({"complete": True, "scored": 5, "elements": [{"id": "E1"}]}))
    assert chain.committable("chr22") is False  # the old shape: megabytes of elements inside
    p.write_text(json.dumps({"complete": False, "scored": 5, "elements_where": "local"}))
    assert chain.committable("chr22") is False  # not finished
    p.write_text(json.dumps({"complete": True, "scored": 5, "elements_where": "local"}))
    assert chain.committable("chr22") is True


def test_scorer_processes_parses_ps_output():
    chain = _chain()
    ps = "\n".join(
        [
            "  PID COMMAND",
            " 100 /usr/bin/Python scripts/enhancer_targets_all.py --chrom chr19",
            " 200 /usr/bin/Python scripts/enhancer_targets_all.py --chrom chr1",
            " 300 /bin/zsh -c pgrep -f 'enhancer_targets_all.py --chrom chr7'",
            " 400 uv run python scripts/enhancer_targets_all_chain.py",
        ]
    )
    assert chain.scorer_processes(ps) == [(100, "chr19"), (200, "chr1")]


def test_suppress_strays(tmp_path, monkeypatch):
    chain = _chain()
    monkeypatch.setattr(chain, "JOBS", tmp_path)
    monkeypatch.setattr(chain, "log", lambda m: None)
    (tmp_path / "enhancer_targets_all_chr1.json").write_text("{}")
    killed = []
    line = " {pid} x Python scripts/enhancer_targets_all.py --chrom {c}"
    ps = line.format(pid=100, c="chr19") + chr(10) + line.format(pid=200, c="chr1")
    out = chain.suppress_strays("chr19", ps, kill=lambda pid, sig: killed.append((pid, sig)))
    assert out == ["chr1"] and [p for p, _ in killed] == [200]
    assert (tmp_path / "enhancer_targets_all_chr1.superseded-by-chain").exists()
    assert not (tmp_path / "enhancer_targets_all_chr1.json").exists()


def _ccres(path, dels, pels, pls=1):
    import gzip

    with gzip.open(path, "wt") as fh:
        fh.write("# header\n")
        for cls, n in (("dELS", dels), ("pELS", pels), ("PLS", pls)):
            for i in range(n):
                fh.write(f"chrT\t{i}\t{i + 10}\tE{cls}{i}\t{cls}\t0\n")


def test_a_capped_run_is_not_complete(tmp_path, monkeypatch):
    chain = _chain()
    monkeypatch.setattr(chain, "RESULTS", tmp_path)
    logged = []
    monkeypatch.setattr(chain, "log", logged.append)
    _ccres(tmp_path / "ccres_chrT.bed.gz", dels=3, pels=2)
    assert chain.expected_elements("chrT") == 5
    result = tmp_path / "enhancer_targets_all_chrT.json"
    result.write_text(json.dumps({"elements_total": 2, "scored": 2, "complete": True}))
    assert chain.complete("chrT") is False and "covers 2 of 5" in logged[-1]
    result.write_text(json.dumps({"elements_total": 5, "scored": 5, "complete": True}))
    assert chain.complete("chrT") is True
    result.write_text(json.dumps({"elements_total": 5, "scored": 4, "complete": True}))
    assert chain.complete("chrT") is False


def test_stall_state_flags_a_frozen_marker():
    chain = _chain()
    last, since, stalled = chain.stall_state(100.0, None, 0.0, 10.0, limit=600)
    assert (last, since, stalled) == (100.0, 10.0, False)
    same = chain.stall_state(100.0, 100.0, 10.0, 300.0, limit=600)
    assert same == (100.0, 10.0, False)  # frozen, but not yet past the limit
    stuck = chain.stall_state(100.0, 100.0, 10.0, 700.0, limit=600)
    assert stuck == (100.0, 10.0, True)  # frozen for 690 s: stalled


def test_progress_marker_reads_heartbeat_and_result(tmp_path, monkeypatch):
    chain = _chain()
    monkeypatch.setattr(chain, "JOBS", tmp_path)
    monkeypatch.setattr(chain, "RESULTS", tmp_path)
    assert chain.progress_marker("chrT") is None
    beat = tmp_path / "enhancer_targets_all_chrT.heartbeat"
    beat.write_text("1")
    assert chain.progress_marker("chrT") == beat.stat().st_mtime
