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
