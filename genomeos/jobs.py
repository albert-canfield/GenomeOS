"""Background jobs started from the UI or CLI, with progress visible.

A job is a subprocess writing a log under data/jobs/<name>.log. Progress is
read from the log and, where a job writes a result file incrementally, from
that file. Only whitelisted jobs can be started (no arbitrary commands).
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

JOBS_DIR = Path("data/jobs")


# name -> (argv, total-steps hint, result name whose keys count progress)
_STEP = re.compile(r"(chr[0-9XYM]+): (\d+)/(\d+) compiled")


def _current_step(name: str, root: Path) -> tuple[str, int, int] | None:
    """The newest 'chrN: a/b compiled' line of a job's log (the live process writes the latest one)."""
    log = root / "data" / "jobs" / f"{name}.log"
    if not log.exists():
        return None
    best = None
    for m in _STEP.finditer(log.read_text(errors="replace")):
        best = (m.group(1), int(m.group(2)), int(m.group(3)))
    return best


def _proteome_progress(root: Path) -> float:
    done = len(list((root / "data" / "results").glob("proteome_chr*.json")))
    step = _current_step("proteome_genome_wide", root)
    if step and not (root / "data" / "results" / f"proteome_{step[0]}.json").exists() and step[2]:
        return round(done + step[1] / step[2], 2)
    return float(done)


def _result_count(root: Path, name: str, key: str) -> int:
    p = root / "data" / "results" / f"{name}.json"
    try:
        return int(json.loads(p.read_text()).get("summary", {}).get(key, 0)) if p.exists() else 0
    except (OSError, json.JSONDecodeError, ValueError):
        return 0


def _count_curated(root: Path) -> int:
    n = 0
    for p in (root / "data" / "results").glob("unknown_chr*.json"):
        try:
            if json.loads(p.read_text()).get("curated_repeats"):
                n += 1
        except (OSError, json.JSONDecodeError):
            continue
    return n


CATALOG: dict[str, dict] = {
    "anatomy_genome_wide": {
        "argv": [sys.executable, "scripts/anatomy_genome_wide.py"],
        "describe": "Every human chromosome: download, count the blocks, delete; keep one summary.",
        "total": 25,
        "result": "anatomy_hg38_by_chromosome",
        "count": lambda r: len(r.get("chromosomes", {})),
    },
    "signals_learn_chr21": {
        "argv": [sys.executable, "-m", "genomeos.cli", "signals", "learn", "--chrom", "chr21"],
        "describe": "Learn the splice and start signals from chromosome 21.",
        "total": 1,
        "result": "signals_chr21",
        "count": lambda r: 1 if r else 0,
    },
    "unknown_chr21": {
        "argv": [sys.executable, "-m", "genomeos.cli", "unknown", "--chrom", "chr21"],
        "describe": "Classify every UNKNOWN block of chromosome 21, largest first.",
        "total": 1,
        "result": "unknown_chr21",
        "count": lambda r: 1 if r else 0,
    },
    "unknown_genome_wide": {
        "argv": [sys.executable, "scripts/unknown_genome_wide.py"],
        "describe": "Every chromosome: classify the UNKNOWN blocks with curated repeats; keep the summary.",
        "total": 24,
        "result": "unknown_genome_wide",
        # progress = chromosomes classified with the curated repeat pass (a re-run starts from the
        # chromosomes that lack it, so the earlier pass does not count)
        "progress": lambda root: _count_curated(root),
        "count": None,
        "complete": lambda root: _count_curated(root) >= 24,
        "auto_heal": True,
    },
    "proteome_chr21": {
        "argv": [sys.executable, "-m", "genomeos.cli", "proteome", "--chrom", "chr21"],
        "describe": "Compile every protein of chromosome 21 from seven public databases; keep the coverage.",
        "total": 1,
        "result": "proteome_chr21",
        "count": lambda r: 1 if r else 0,
    },
    "proteome_genome_wide": {
        "argv": [sys.executable, "scripts/proteome_genome_wide.py"],
        "describe": "Every chromosome's proteins compiled from seven public databases, smallest first.",
        "total": 25,
        "result": None,
        "count": None,
        "progress": lambda root: _proteome_progress(root),
        "complete": lambda root: len(list((root / "data" / "results").glob("proteome_chr*.json"))) >= 25,
        "auto_heal": True,
    },
    "distil": {
        "argv": [sys.executable, "-m", "genomeos.cli", "data", "distil"],
        "describe": "Turn any raw downloads present into result summaries.",
        "total": 1,
        "result": None,
        "count": None,
    },
}

CATALOG["fetch_hg002"] = {
    "argv": [sys.executable, "-m", "genomeos.cli", "data", "fetch", "--individual", "--chrom", "chr22"],
    "describe": "Stream the GIAB HG002 benchmark once; keep each chromosome's variants for twins/lookups.",
    "total": 1,
    "result": None,
    "count": None,
}

# one fetch job per human chromosome: bring it to full footing (sequence, models, elements, repeats)
CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]


def _enhancer_job(chrom: str, sample: int = 200) -> dict:
    """AlphaGenome feature b on one chromosome: a sampled, resumable, self-healing job (needs the key)."""
    name = f"enhancer_targets_{chrom}"
    return {
        "argv": [sys.executable, "scripts/enhancer_targets.py", "--chrom", chrom, "--sample", str(sample)],
        "describe": f"AlphaGenome: delete {sample} {chrom} distal enhancers one at a time; which gene moves?",
        "total": sample,
        "result": name,
        "count": lambda r: r.get("summary", {}).get("elements_scored", 0),
        "complete": lambda root, n=name, s=sample: _result_count(root, n, "elements_scored") >= s,
        "auto_heal": True,
    }


for _c in CHROMOSOMES:
    if _c != "chrM":
        CATALOG[f"enhancer_targets_{_c}"] = _enhancer_job(_c)
for _c in CHROMOSOMES:
    CATALOG[f"fetch_{_c}"] = {
        "argv": [sys.executable, "-m", "genomeos.cli", "data", "fetch", "--analyse", "--chrom", _c],
        "describe": f"Fetch {_c}: sequence, GENCODE rows, ENCODE elements, RepeatMasker.",
        "total": 1,
        "result": f"rmsk_{_c}",
        "count": lambda r: 1 if r else 0,
    }


@dataclass(slots=True)
class JobStatus:
    name: str
    describe: str
    state: str  # idle | running | done | failed
    started: float | None
    finished: float | None
    done: float
    total: int
    last_lines: list[str]
    detail: str = ""
    activity: float | None = None  # seconds since the last sign of life
    heals: int = 0

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "describe": self.describe,
            "state": self.state,
            "started": self.started,
            "finished": self.finished,
            "done": self.done,
            "total": self.total,
            "fraction": round(self.done / self.total, 3) if self.total else None,
            "last_lines": self.last_lines,
            "detail": self.detail,
            "activity": self.activity,
            "heals": self.heals,
        }


_running: dict[str, subprocess.Popen] = {}


def _meta_path(name: str) -> Path:
    return JOBS_DIR / f"{name}.json"


def _alive(pid: int | None) -> bool:
    """Is a job process started by any earlier server or shell still running?"""
    if not pid:
        return False
    try:
        os.kill(pid, 0)
    except (ProcessLookupError, PermissionError):
        return False
    return True


STALL_AFTER = 20 * 60  # seconds without a heartbeat or a log line before a running job counts as stalled
HEAL_COOLDOWN = 10 * 60  # seconds between automatic restarts of the same job


def heartbeat(name: str, root: Path = Path(".")) -> None:
    """Called by a job script at every step: proof of life for the supervisor."""
    p = root / "data" / "jobs" / f"{name}.heartbeat"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(str(time.time()))


def last_activity(name: str, root: Path = Path(".")) -> float | None:
    """Seconds since the job last showed life (heartbeat file or log write), None if never."""
    stamps = []
    for p in (root / "data" / "jobs" / f"{name}.heartbeat", root / "data" / "jobs" / f"{name}.log"):
        if p.exists():
            stamps.append(p.stat().st_mtime)
    return round(time.time() - max(stamps)) if stamps else None


def is_complete(name: str, root: Path = Path(".")) -> bool:
    fn = CATALOG[name].get("complete")
    return bool(fn(root)) if fn else False


def heal(name: str, root: Path = Path("."), force: bool = False) -> JobStatus:
    """Bring a job back on track: if it is not complete and its process is gone (or stalled), start it
    again; the scripts resume from what is already saved. A live, active job is left alone."""
    st = status(name, root)
    if is_complete(name, root) and not force:
        return st
    pid = _recorded_pid(name)
    if st.state == "running" and not force:
        return st
    if st.state == "stalled" and _alive(pid):
        try:
            os.kill(pid, 15)
            time.sleep(2)
        except ProcessLookupError:
            pass
    _running.pop(name, None)
    return start(name, root)


def supervise_once(root: Path = Path(".")) -> list[str]:
    """One supervisor pass: restart every auto-heal job that was started, is not complete, and is
    dead or stalled, at most once per HEAL_COOLDOWN. Returns the names healed."""
    healed = []
    for name, spec in CATALOG.items():
        if not spec.get("auto_heal") or not _meta_path(name).exists() or is_complete(name, root):
            continue
        st = status(name, root)
        if st.state not in ("failed", "unknown", "stalled"):
            continue
        meta = json.loads(_meta_path(name).read_text())
        if time.time() - (meta.get("healed") or 0) < HEAL_COOLDOWN:
            continue
        heal(name, root)
        meta = json.loads(_meta_path(name).read_text())
        meta["healed"] = time.time()
        meta["heals"] = (meta.get("heals") or 0) + 1
        _meta_path(name).write_text(json.dumps(meta))
        healed.append(name)
    return healed


def supervise(root: Path = Path("."), interval: int = 60) -> None:
    """Run forever (a daemon thread in the server): heal what needs healing every `interval` seconds."""
    import contextlib

    while True:
        with contextlib.suppress(Exception):  # the supervisor itself must never die
            supervise_once(root)
        time.sleep(interval)


def _recorded_pid(name: str) -> int | None:
    if not _meta_path(name).exists():
        return None
    meta = json.loads(_meta_path(name).read_text())
    return meta.get("pid") if meta.get("finished") is None else None


def start(name: str, root: Path = Path(".")) -> JobStatus:
    if name not in CATALOG:
        raise KeyError(f"unknown job {name!r}; known: {sorted(CATALOG)}")
    if name in _running and _running[name].poll() is None:
        return status(name, root)
    if _alive(_recorded_pid(name)):
        # started by an earlier server process (the server restarts; the job does not): do not start a twin
        return status(name, root)
    JOBS_DIR.mkdir(parents=True, exist_ok=True)
    log = JOBS_DIR / f"{name}.log"
    log.write_text("")
    # append mode: every write lands at the end, so a log can never carry stale lines from an
    # earlier process after its own newer ones
    fh = open(log, "a")  # noqa: SIM115  (the subprocess owns the handle)
    proc = subprocess.Popen(CATALOG[name]["argv"], cwd=root, stdout=fh, stderr=subprocess.STDOUT)  # noqa: S603
    _running[name] = proc
    _meta_path(name).write_text(
        json.dumps({"pid": proc.pid, "started": time.time(), "finished": None, "code": None})
    )
    return status(name, root)


def status(name: str, root: Path = Path(".")) -> JobStatus:
    spec = CATALOG[name]
    meta = json.loads(_meta_path(name).read_text()) if _meta_path(name).exists() else {}
    proc = _running.get(name)
    state = "idle"
    if proc is not None:
        code = proc.poll()
        if code is None:
            state = "running"
        else:
            state = "done" if code == 0 else "failed"
            if meta.get("finished") is None:
                meta.update({"finished": time.time(), "code": code})
                _meta_path(name).write_text(json.dumps(meta))
    elif meta:
        if meta.get("code") is None and _alive(meta.get("pid")):
            state = "running"
        else:
            state = (
                "done" if meta.get("code") == 0 else "failed" if meta.get("code") is not None else "unknown"
            )
    activity = last_activity(name, root)
    if state == "running" and activity is not None and activity > STALL_AFTER:
        state = "stalled"
    if state in ("unknown", "failed") and is_complete(name, root):
        state = "done"  # the results say it finished, whatever happened to the process that made them
    done: float = 0
    detail = ""
    if spec.get("progress"):
        done = spec["progress"](root)
        step = _current_step(name, root)
        if step:
            detail = f"{step[0]}: {step[1]:,}/{step[2]:,}"
    if spec.get("result") and spec.get("count"):
        from genomeos.results import load_result

        rp = root / "data" / "results" / f"{spec['result']}.json"
        r = load_result(spec["result"], root / "data" / "results")
        done = spec["count"](r) if r else 0
        # a job started outside this process (CLI, nohup) shows as running while its result keeps changing
        if (
            state in ("idle", "unknown")
            and rp.exists()
            and done < spec["total"]
            and time.time() - rp.stat().st_mtime < 600
        ):
            state = "running"
        if state in ("idle", "unknown") and done >= spec["total"]:
            state = "done"
    elif state == "done":
        done = spec["total"]
    log = JOBS_DIR / f"{name}.log"
    lines = log.read_text().splitlines()[-6:] if log.exists() else []
    return JobStatus(
        name,
        spec["describe"],
        state,
        meta.get("started"),
        meta.get("finished"),
        done,
        spec["total"],
        lines,
        detail=detail,
        activity=activity,
        heals=int(meta.get("heals") or 0) if meta else 0,
    )


def all_status(root: Path = Path(".")) -> list[JobStatus]:
    return [status(name, root) for name in CATALOG]
