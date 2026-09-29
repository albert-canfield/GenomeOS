"""The two C. elegans writers give the same output under different string-hash seeds (lane-repro, 2026-09-29).

`scripts/celegans_commitment.py` did not: `competence_test` iterated a set of factor names, so the order of
its rows, hence of its strata, hence which seeded shuffle the null drew for which stratum, followed
PYTHONHASHSEED. Between seeds 0 and 1, 12 values of the competence block moved in the fourth decimal and
the committed file matched neither. The fix sorts the set; before it, the first test below failed on the
committed code: 12 values differed at the writer's full scale and 14 at `--permutations 2`.

Each writer runs in a subprocess per seed, with `save_result` captured so nothing is written into the
checkout, and with PYTHONPATH set to this checkout: in a separate git worktree the editable install would
otherwise import the other checkout's `genomeos`. The fate-reads writer runs its harvest, the threshold
plateau and one arm (the shipped read, in sample, under both recheck modes) rather than all 17 arms and 88
folds, which take five minutes; the full runs were compared by hand when the results were regenerated.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SEEDS = ("0", "1")
NEEDS = [
    ROOT / "data/results/celegans_lineage_cells.json",
    ROOT / "data/results/celegans_tf_atlas_cells.json",
    ROOT / "data/knowledge/celegans/atlas_levels.json",
]

COMMITMENT = """
import contextlib, importlib.util, io, json, sys
spec = importlib.util.spec_from_file_location("writer", "scripts/celegans_commitment.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
out = {}
m.save_result = lambda name, payload, **kw: out.update(payload)
sys.argv = ["celegans_commitment.py", "--permutations", "2"]
with contextlib.redirect_stdout(io.StringIO()):
    m.main()
out.pop("date", None)
print(json.dumps(out, sort_keys=True))
"""

FATE_READS = """
import importlib.util, json, shutil, tempfile
from pathlib import Path
spec = importlib.util.spec_from_file_location("writer", "scripts/celegans_fate_reads.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
ref = m.ReferenceLineage.load()
labels = {c.id: c.tissue for c in m.fr.embryonic_terminal(ref)}
cells = sorted(labels)
factors = sorted({f for fs in m.load_cells().values() for f in fs})
reads, points = m.harvest(ref, factors)
out = {"points": {str(k): v for k, v in points.items()}, "reads": reads}
out["plateau"] = m.plateau(reads, labels, m.presence_at_birth(ref))
feats = m.fr.runtime_features(reads, "exposure", "lineage", m.fr.THRESHOLD_MIN)
out["fates.bio"] = m.fr.to_bio_fates(m.program_rules(feats, labels, cells), read="exposure", window="lineage")
tmp = Path(tempfile.mkdtemp(prefix="determinism-"))
try:
    arm = m.Arm(tmp / "arm")

    def write(fit):
        rules = m.program_rules(feats, labels, fit)
        (arm.dir / "fates.bio").write_text(m.fr.to_bio_fates(rules, read="exposure", window="lineage"))

    arm_out = m.score_arm(arm, ref, write, labels, folds=False)
    arm_out["in_sample"].pop("to_adult")  # the run to 6,000 min is the slow part and reads nothing new
    out["arm"] = arm_out
finally:
    shutil.rmtree(tmp, ignore_errors=True)
print(json.dumps(out, sort_keys=True))
"""


def _run_under_seeds(code: str) -> dict[str, str]:
    """Start one subprocess per hash seed, together, and return each one's stdout."""
    procs = {}
    for seed in SEEDS:
        env = {**os.environ, "PYTHONHASHSEED": seed, "PYTHONPATH": str(ROOT)}
        procs[seed] = subprocess.Popen(
            [sys.executable, "-c", code],
            cwd=ROOT,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    out = {}
    for seed, p in procs.items():
        stdout, stderr = p.communicate(timeout=600)
        assert p.returncode == 0, f"seed {seed}: {stderr[-2000:]}"
        out[seed] = stdout.strip().splitlines()[-1]
    return out


@pytest.mark.skipif(not all(p.exists() for p in NEEDS), reason="the worm's distilled inputs are not here")
def test_commitment_writer_does_not_depend_on_the_hash_seed():
    runs = _run_under_seeds(COMMITMENT)
    a, b = (json.loads(runs[s]) for s in SEEDS)
    assert a["competence"]["all_factors"]["onsets"] == 15353  # the real data, not an empty run
    assert runs[SEEDS[0]] == runs[SEEDS[1]]


@pytest.mark.skipif(not all(p.exists() for p in NEEDS), reason="the worm's distilled inputs are not here")
def test_fate_reads_writer_does_not_depend_on_the_hash_seed():
    runs = _run_under_seeds(FATE_READS)
    a = json.loads(runs[SEEDS[0]])
    assert len(a["reads"]) == 555 and a["arm"]["in_sample"]["fates_checked"] == 555
    assert runs[SEEDS[0]] == runs[SEEDS[1]]
