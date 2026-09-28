# SPDX-License-Identifier: AGPL-3.0-or-later
"""No code writes into data/results/ except through `save_result` (item 12 S6 follow-up, lane-contract).

`save_result` (genomeos/results.py) is the result contract: a new result without a complete manifest is
quarantined, never written. A writer that opens a file in data/results/ itself is outside it, which is
how 20 writers came to put or edit results there unchecked (census d6a77dd). This test reads every
Python file under genomeos/ and scripts/ with scripts/results_writer_census.py's static analysis and
fails on any site that writes into data/results/ other than through save_result, unless it is on one of
two short lists below, each entry with its reason. An entry that no longer matches a site fails too, so
the lists cannot outlive what they excuse.

The analysis is checked on its own first: every pattern the census found in the tree, written small,
must be caught, and the reads and the sanctioned write must not.
"""

from __future__ import annotations

import importlib.util
import sys
import textwrap
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
_SPEC = importlib.util.spec_from_file_location(
    "results_writer_census", ROOT / "scripts/results_writer_census.py"
)
census = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = census  # its dataclass resolves annotations through sys.modules
_SPEC.loader.exec_module(census)

#: Writes into a data/results/ directory that are not a write into this registry. (file, function): why.
EXCEPTIONS = {
    ("scripts/manifest_rebuild.py", "rebuild"): "copies a rebuilt result under another name inside the "
    "disposable second checkout (<where>/rebuild-<sha>/data/results), never this registry",
    ("scripts/activator_combination.py", "main"): "writes only when --out is given, to args.out; the "
    "registry write is save_result (a false positive of the flow-insensitive analysis)",
}

#: Fetched reference data kept in data/results/ by an older layout. Not results: none is JSON, the registry
#: lists *.json only, and save_result's manifest is a key inside a result's JSON. Out of scope for the
#: result contract; moving them out of data/results is a separate change (docs/DATA.md).
OUT_OF_SCOPE_DATA = {
    ("genomeos/genome/duplications.py", "save_rows"): "superdups_<chrom>.bed.gz",
    ("genomeos/genome/fetch.py", "fetch_gencode_chrom"): "gencode_v50_<chrom>.gff3.gz",
    ("genomeos/genome/mouse.py", "fetch_ccres"): "ccres_mm10_<chrom>.bed.gz",
    ("genomeos/genome/mouse.py", "fetch_orthology"): "mgi_mouse_human_orthology.tsv.gz",
    ("genomeos/genome/mouse.py", "fetch_compara_orthology"): "compara_mouse_human_orthology.tsv.gz",
    ("genomeos/genome/reader.py", "fetch_peaks"): "dnase_<biosample>_<chrom>.bed.gz",
    ("genomeos/genome/regulatory.py", "save_ccres"): "ccres_<chrom>.bed.gz",
    ("genomeos/genome/repeats.py", "save_repeats"): "rmsk_<chrom>.bed.gz",
    ("genomeos/storage.py", "_distil_hg002_chr21"): "HG002_chr21.vcf",
    ("genomeos/storage.py", "_distil_gencode_subset"): "gencode_v50_chr21_chrM.gff3.gz",
    ("scripts/unknown_genome_wide.py", "main"): "deletes rmsk_<chrom>.bed.gz",
}


def _flagged(tmp_path: Path, files: dict[str, str]) -> list:
    for rel, text in files.items():
        p = tmp_path / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(textwrap.dedent(text))
    return census.World(tmp_path, list(files)).sites()


POSITIVE = {
    "a string path": """
        from pathlib import Path
        Path("data/results/x.json").write_text("{}")
    """,
    "a path joined from data and results": """
        from pathlib import Path
        ROOT = Path(__file__).resolve().parents[1]
        OUT = ROOT / "data" / "results" / "x.json"
        def main():
            OUT.write_text("{}")
    """,
    "os.path.join and open for writing": """
        import json, os
        RESULTS = os.path.join(os.path.dirname(__file__), "data", "results")
        def main():
            with open(os.path.join(RESULTS, "x.json"), "w") as f:
                json.dump({}, f)
    """,
    "the registry constant": """
        from genomeos.results import RESULTS_DIR
        def main():
            (RESULTS_DIR / "x.json").write_text("{}")
    """,
    "a writer helper": """
        from genomeos.results import RESULTS_DIR
        def _write(path, text):
            path.write_text(text)
        def main():
            _write(RESULTS_DIR / "x.json", "{}")
    """,
    "a default argument": """
        from pathlib import Path
        def save(table, path=Path("data/results/t.json")):
            path.write_text(str(table))
    """,
    "an argparse default": """
        import argparse
        from pathlib import Path
        def main():
            ap = argparse.ArgumentParser()
            ap.add_argument("--out", default="data/results/x.json")
            args = ap.parse_args()
            Path(args.out).write_text("{}")
    """,
    "a copy into the registry": """
        import shutil
        from genomeos.results import RESULTS_DIR
        def main(src):
            shutil.copyfile(src, RESULTS_DIR / "x.json")
    """,
    "a rename into the registry": """
        import os
        from genomeos.results import RESULTS_DIR
        def main(tmp):
            os.replace(tmp, RESULTS_DIR / "x.json")
    """,
    "an in-place edit over a glob": """
        import json
        from genomeos.results import RESULTS_DIR
        def main():
            for p in RESULTS_DIR.glob("*.json"):
                d = json.loads(p.read_text())
                d["note"] = 1
                p.write_text(json.dumps(d))
    """,
    "gzip into the registry": """
        import gzip
        from genomeos.results import RESULTS_DIR
        def main():
            with gzip.open(RESULTS_DIR / "x.json.gz", "wt") as f:
                f.write("{}")
    """,
    "a method on an object": """
        from pathlib import Path
        class Store:
            def save(self, path):
                Path(path).write_text("{}")
        def main():
            Store().save("data/results/s.json")
    """,
    "a results_dir parameter": """
        def run(results_dir):
            (results_dir / "x.json").write_text("{}")
    """,
}


@pytest.mark.parametrize("name", sorted(POSITIVE))
def test_the_analysis_catches_each_pattern_the_census_found(tmp_path, name):
    sites = _flagged(tmp_path, {"scripts/w.py": POSITIVE[name]})
    assert [s for s in sites if s.kind in ("write", "possible")], name


def test_the_analysis_follows_an_imported_constant_and_a_path_function(tmp_path):
    sites = _flagged(
        tmp_path,
        {
            "genomeos/a.py": 'from pathlib import Path\nCELLS = Path("data/results/c.json")\n'
            'def out(name):\n    return Path("data/results") / f"{name}.json"\n',
            "scripts/b.py": "from genomeos.a import CELLS, out\ndef main():\n    CELLS.write_text('')\n"
            "    out('y').write_text('')\n",
        },
    )
    assert sorted(s.line for s in sites if s.file == "scripts/b.py") == [3, 4]


NEGATIVE = """
    import json, subprocess
    from pathlib import Path
    from genomeos.results import RESULTS_DIR, save_result
    def main(m):
        save_result("x", {"a": 1}, manifest=m)
        Path("data/cache/x.json").write_text("{}")
        Path("data/results_legacy.txt").write_text("")
        d = json.loads((RESULTS_DIR / "x.json").read_text())
        with open(RESULTS_DIR / "x.json") as f:
            d = json.load(f)
        subprocess.run(["git", "show", "HEAD:data/results/x.json"], check=True)
        return d
"""


def test_reads_the_sanctioned_write_and_other_directories_are_not_flagged(tmp_path):
    assert _flagged(tmp_path, {"scripts/n.py": NEGATIVE}) == []


def test_nothing_under_genomeos_or_scripts_writes_into_data_results_but_save_result():
    sites = census.find_sites(ROOT)
    allowed = set(EXCEPTIONS) | set(OUT_OF_SCOPE_DATA)
    outside = [
        f"{s.file}:{s.line} {s.function} {s.sink}({s.target}) {s.via}".strip()
        for s in sites
        if (s.file, s.function) not in allowed
    ]
    assert outside == [], "write into data/results through save_result, or justify it here:\n" + "\n".join(
        outside
    )


def test_every_listed_exception_still_matches_a_site():
    sites = {(s.file, s.function) for s in census.find_sites(ROOT)}
    stale = sorted(k for k in (*EXCEPTIONS, *OUT_OF_SCOPE_DATA) if k not in sites)
    assert stale == [], f"remove these entries, they no longer excuse anything: {stale}"
