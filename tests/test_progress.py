from pathlib import Path

from genomeos import jobs
from genomeos.web.server import Api


def test_progress_rows_come_from_the_docs():
    p = Api(Path(".")).progress()
    assert p["count"] > 20 and p["done"] > 10
    assert all({"task", "status", "evidence"} <= set(r) for r in p["rows"])


def test_job_catalog_and_status():
    names = {j.name for j in jobs.all_status(Path("."))}
    assert {"anatomy_genome_wide", "signals_learn_chr21", "distil"} <= names
    st = jobs.status("anatomy_genome_wide", Path("."))
    assert (
        st.total == 25 and 0 <= st.done <= 25 and st.state in ("idle", "running", "done", "failed", "unknown")
    )
    d = st.to_dict()
    assert "fraction" in d and "last_lines" in d


def test_unknown_job_is_rejected():
    import pytest

    with pytest.raises(KeyError):
        jobs.start("rm_rf", Path("."))


def test_all_elements_job_reads_its_result_not_its_process(tmp_path, monkeypatch):
    """A chromosome the chain has finished reads done, with the element count its result records.

    The scorer is stopped by the chain rather than exiting on its own, so the job's metadata says
    "started, never finished" for a run that did complete. Before this, finished chromosomes showed
    as stalled and their totals came from whatever the result held when the module was imported.
    """
    import json
    from pathlib import Path

    from genomeos import jobs

    results = tmp_path / "data" / "results"
    results.mkdir(parents=True)
    (results / "enhancer_targets_all_chr21.json").write_text(
        json.dumps({"chrom": "chr21", "scored": 12139, "elements_total": 12139, "complete": True})
    )
    (results / "enhancer_targets_all_chr15.json").write_text(
        json.dumps({"chrom": "chr15", "scored": 3400, "elements_total": 31058, "complete": False})
    )
    meta = tmp_path / "data" / "jobs"
    meta.mkdir(parents=True)
    monkeypatch.setattr(jobs, "JOBS_DIR", meta)
    # a process the chain stopped mid-run: a pid far above any in use, so it reads as gone
    for name in ("enhancer_targets_all_chr21", "enhancer_targets_all_chr15"):
        (meta / f"{name}.json").write_text(json.dumps({"pid": 2**22, "started": 1.0, "finished": None}))

    done = jobs.status("enhancer_targets_all_chr21", Path(tmp_path))
    assert done.state == "done" and done.done == 12139 and done.total == 12139
    part = jobs.status("enhancer_targets_all_chr15", Path(tmp_path))
    assert part.state != "done" and part.done == 3400 and part.total == 31058  # total read per call
    assert jobs.CATALOG["enhancer_targets_all_chr15"]["auto_heal"] is False  # the chain sequences these
