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
from typing import Any

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

CATALOG["budget_genome_wide"] = {
    "argv": [sys.executable, "scripts/budget_genome_wide.py"],
    "describe": "The 98%: Zoonomia constraint over every UNKNOWN block and a best guess per block.",
    "total": 24,
    "result": None,
    "count": None,
    "progress": lambda root: len(list((root / "data" / "results").glob("budget_chr*.json"))),
    "complete": lambda root: len(list((root / "data" / "results").glob("budget_chr*.json"))) >= 24,
    "auto_heal": True,
}

CATALOG["satmut_vista"] = {
    "argv": [sys.executable, "scripts/satmut_vista.py", "--run", "--quota-handed"],
    "describe": "VISTA in-silico mutagenesis: 3,200 windows over the panel's two arms (pre-registered).",
    "total": 3200,
    "result": "satmut_vista",
    "count": lambda r: r.get("windows_scored", 0) if r else 0,
}

CATALOG["executor_e1_extension"] = {
    "argv": [
        sys.executable,
        "scripts/executor_test.py",
        "--wide",
        "--run",
        "--quota-handed",
        "--extension",
        "--max-pairs",
        "2260",
    ],
    "describe": "E1's remaining 2,260 MPRA allele pairs, their own sample (pre-registered 2026-09-16).",
    "total": 2260,
    "result": "executor_mpra_wide_extension",
    "count": lambda r: r.get("pairs_done", 0) if r else 0,
}

CATALOG["human_panel_sweep"] = {
    "argv": [sys.executable, "scripts/human_panel_sweep.py"],
    "describe": "The HPRC panel of 90 human genomes, one unread chromosome at a time (no model key).",
    "total": 24,
    "result": None,
    "count": None,
    "progress": lambda root: len(list((root / "data" / "results").glob("human_panel_chr*.json"))),
    "complete": lambda root: len(list((root / "data" / "results").glob("human_panel_chr*.json"))) >= 24,
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


def _enhancer_chromosomes_done(root: Path, sample: int = 200) -> int:
    return sum(
        1
        for c in CHROMOSOMES
        if c != "chrM" and _result_count(root, f"enhancer_targets_{c}", "elements_scored") >= sample
    )


CATALOG["enhancer_targets_genome_wide"] = {
    "argv": [sys.executable, "scripts/enhancer_targets_genome_wide.py"],
    "describe": "AlphaGenome: 200 distal enhancers per chromosome, one chromosome at a time, until all 24.",
    "total": 24,
    "result": None,
    "count": None,
    "progress": lambda root: float(_enhancer_chromosomes_done(root)),
    "complete": lambda root: _enhancer_chromosomes_done(root) >= 24,
    "auto_heal": True,
}


def _constrained_chromosomes_done(root: Path) -> int:
    return sum(
        1
        for c in CHROMOSOMES
        if c != "chrM" and _result_count(root, f"constrained_targets_{c}", "elements_scored") > 0
    )


CATALOG["constrained_targets_genome_wide"] = {
    "argv": [sys.executable, "scripts/constrained_targets_genome_wide.py", "--run"],
    "describe": "AlphaGenome on the 100 most constrained distal enhancers per chromosome, then one summary.",
    "total": 24,
    "result": "constrained_targets_genome_wide",
    "count": None,
    "progress": lambda root: float(_constrained_chromosomes_done(root)),
    "complete": lambda root: _constrained_chromosomes_done(root) >= 24,
    "auto_heal": True,
}


def _vista_chromosomes_done(root: Path) -> int:
    return sum(1 for c in CHROMOSOMES if (root / "data" / "results" / f"vista_{c}.json").exists())


CATALOG["vista_genome_wide"] = {
    "argv": [sys.executable, "scripts/vista_genome_wide.py", "--run"],
    "describe": "VISTA's measured enhancers against registry, node, constraint and deletion, per chromosome.",
    "total": 23,
    "result": "vista_genome_wide",
    "count": None,
    "progress": lambda root: float(_vista_chromosomes_done(root)),
    "complete": lambda root: _vista_chromosomes_done(root) >= 23,
    "auto_heal": True,
}


def _mpra_chromosomes_done(root: Path) -> int:
    return sum(1 for c in CHROMOSOMES if (root / "data" / "results" / f"mpra_{c}.json").exists())


CATALOG["mpra_genome_wide"] = {
    "argv": [sys.executable, "scripts/mpra_genome_wide.py", "--run"],
    "describe": "lentiMPRA activity in K562, HepG2 and WTC11 against registry, reader, constraint and model.",
    "total": 24,
    "result": "mpra_genome_wide",
    "count": None,
    "progress": lambda root: float(_mpra_chromosomes_done(root)),
    "complete": lambda root: _mpra_chromosomes_done(root) >= 24,
    "auto_heal": True,
}


def _all_elements_total(root: Path, chrom: str) -> int | None:
    """The chromosome's element count once its all-elements result exists (the script writes it first)."""
    p = root / "data" / "results" / f"enhancer_targets_all_{chrom}.json"
    try:
        return int(json.loads(p.read_text()).get("elements_total")) if p.exists() else None
    except (OSError, json.JSONDecodeError, ValueError, TypeError):
        return None


def _all_elements_complete(root: Path, chrom: str) -> bool:
    """The chromosome's own result says the run finished; the chain and the Progress tab both trust it."""
    p = root / "data" / "results" / f"enhancer_targets_all_{chrom}.json"
    try:
        return bool(json.loads(p.read_text()).get("complete")) if p.exists() else False
    except (OSError, json.JSONDecodeError, ValueError, TypeError):
        return False


for _c in CHROMOSOMES:
    if _c == "chrM":
        continue
    CATALOG[f"enhancer_targets_all_{_c}"] = {
        "argv": [sys.executable, "scripts/enhancer_targets_all.py", "--chrom", _c, "--workers", "8"],
        "describe": f"AlphaGenome: every enhancer in a {_c} node deleted, effect per cell; 8 at a time.",
        "total": _all_elements_total(Path("."), _c) or (12139 if _c == "chr21" else None),
        "total_fn": lambda root, _c=_c: _all_elements_total(root, _c),
        "result": f"enhancer_targets_all_{_c}",
        "count": lambda r: int(r.get("scored", 0)),
        "complete": lambda root, _c=_c: _all_elements_complete(root, _c),
        # not auto-healed: each chromosome spends the shared quota and their order matters, so restarts
        # belong to whoever is sequencing them (scripts/enhancer_targets_all_chain.py), not to a
        # supervisor tick that would start a chromosome nobody asked for on the next `serve` restart
        "auto_heal": False,
    }
CATALOG["consequence_targets"] = {
    "argv": [sys.executable, "scripts/consequence_targets.py"],
    "describe": "GWAS lead variants against a shifted control; ClinVar non-coding variants on the elements.",
    "total": None,
    "result": "consequence_targets",
    "count": None,
    "auto_heal": True,
}
CATALOG["eqtl_targets"] = {
    "argv": [sys.executable, "scripts/eqtl_targets.py"],
    "describe": "GTEx eQTLs streamed once (1.4 GB, only the hits kept); both target callers judged.",
    "total": None,
    "result": "eqtl_targets",
    "count": None,
    "auto_heal": True,
}
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


def _is_zombie(pid: int) -> bool:
    """A killed child whose parent has not reaped it: still answers signal 0, runs nothing."""
    try:
        out = subprocess.run(  # noqa: S603, S607 - ps is the portable way to read a process state
            ["ps", "-o", "state=", "-p", str(pid)], capture_output=True, text=True, timeout=5, check=False
        )
    except (OSError, subprocess.SubprocessError):
        return False
    return out.stdout.strip().startswith("Z")


def _alive(pid: int | None) -> bool:
    """Is a job process started by any earlier server or shell still running? A zombie is not: a job
    started through the registry by a driver that never waits leaves one behind, and counting it as
    running refuses every restart of that chromosome (chr19, 2026-09-13)."""
    if not pid:
        return False
    try:
        os.kill(pid, 0)
    except (ProcessLookupError, PermissionError):
        return False
    return not _is_zombie(pid)


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


def adopt(name: str, pid: int, root: Path = Path(".")) -> JobStatus:
    """Record a job that is already running, started outside the registry (a shell, nohup, a driver).

    The Progress tab shows what the registry has a record of, so a job started by hand is invisible
    while it runs and its work looks like nothing happening — which is exactly what it looked like on
    2026-09-16, with the panel reading chr2 and the tab showing an idle lane. Adopting writes the same
    record `start` writes, so the job is shown, stall detection applies to it, and nothing starts a
    twin of it. It refuses a pid that is not running, and one the registry already has alive under
    this name.
    """
    if name not in CATALOG:
        raise KeyError(f"unknown job {name!r}; known: {sorted(CATALOG)}")
    if not _alive(pid):
        raise RuntimeError(f"pid {pid} is not running, so there is nothing to adopt")
    current = _recorded_pid(name)
    if current is not None and current != pid and _alive(current):
        raise RuntimeError(f"{name} already has a live process, pid {current}: two would race")
    JOBS_DIR.mkdir(parents=True, exist_ok=True)
    _meta_path(name).write_text(
        json.dumps({"pid": pid, "started": time.time(), "finished": None, "code": None, "adopted": True})
    )
    return status(name, root)


def _recorded_pid(name: str) -> int | None:
    if not _meta_path(name).exists():
        return None
    meta = json.loads(_meta_path(name).read_text())
    return meta.get("pid") if meta.get("finished") is None else None


# every all-elements job spends the same AlphaGenome key, so two of them only race: a stray started
# from the Progress tab, a supervisor tick or a hand has been killed by the chain four times in a day,
# and each one costs the requests it spends before it dies. The registry refuses the second one instead.
SHARED_KEY_PREFIX = "enhancer_targets_all_"


KEY_LOCK = "alphagenome.key"


def _lock_path() -> Path:
    return JOBS_DIR / f"{KEY_LOCK}.lock"


def key_holder() -> dict[str, Any] | None:
    """Whoever holds the shared model key, or None. A lock whose process is gone is stale and ignored."""
    try:
        held = json.loads(_lock_path().read_text())
    except (OSError, json.JSONDecodeError):
        return None
    if not _alive(held.get("pid")):
        return None
    return held


def take_key(holder: str, what: str = "", pid: int | None = None) -> dict[str, Any]:
    """Claim the shared model key for a job that is not in the registry (a scoring script, a test run).

    Raises if someone else holds it. Release it with `drop_key(holder)`; a holder whose process dies
    releases it by itself, since a lock is only honoured while its pid is alive.
    """
    current = key_holder()
    if current and current.get("holder") != holder:
        raise RuntimeError(
            f"the AlphaGenome key is held by {current.get('holder')} ({current.get('what', '')}), "
            f"started {time.strftime('%H:%M', time.localtime(current.get('started', 0)))}: "
            "one scorer per key, so wait for it or ask whoever is coordinating"
        )
    JOBS_DIR.mkdir(parents=True, exist_ok=True)
    held = {"holder": holder, "what": what, "pid": pid or os.getpid(), "started": time.time()}
    _lock_path().write_text(json.dumps(held))
    return held


def drop_key(holder: str) -> None:
    """Release the shared model key if this holder has it; silent if it does not."""
    current = key_holder()
    if current is None or current.get("holder") == holder:
        _lock_path().unlink(missing_ok=True)


def _key_held_by(name: str) -> str | None:
    """Whoever holds the shared model key right now, if this job would contend for it.

    Two things can hold it: another all-elements job in the registry, or a job outside the registry
    that took the lock (the executor runs, a panel scoring pass). Before the lock, a Start from the
    Progress tab during an executor run shared the quota silently for six minutes.
    """
    if not name.startswith(SHARED_KEY_PREFIX):
        return None
    held = key_holder()
    if held and held.get("holder") != name:
        return f"{held.get('holder')} ({held.get('what', '')})"
    for other in CATALOG:
        if other == name or not other.startswith(SHARED_KEY_PREFIX):
            continue
        if (other in _running and _running[other].poll() is None) or _alive(_recorded_pid(other)):
            return other
    return None


def start(name: str, root: Path = Path(".")) -> JobStatus:
    if name not in CATALOG:
        raise KeyError(f"unknown job {name!r}; known: {sorted(CATALOG)}")
    if name in _running and _running[name].poll() is None:
        return status(name, root)
    if _alive(_recorded_pid(name)):
        # started by an earlier server process (the server restarts; the job does not): do not start a twin
        return status(name, root)
    held = _key_held_by(name)
    if held:
        how = (
            "the chain (scripts/enhancer_targets_all_chain.py) sequences the chromosomes one at a "
            "time, so start them through it rather than singly"
            if held.startswith(SHARED_KEY_PREFIX)
            else "wait for it to finish, or ask whoever is coordinating to hand the key over"
        )
        raise RuntimeError(f"{name} shares one AlphaGenome key with {held}, which is running: {how}")
    JOBS_DIR.mkdir(parents=True, exist_ok=True)
    log = JOBS_DIR / f"{name}.log"
    log.write_text("")
    # append mode: every write lands at the end, so a log can never carry stale lines from an
    # earlier process after its own newer ones
    fh = open(log, "a")  # noqa: SIM115  (the subprocess owns the handle)
    proc = subprocess.Popen(CATALOG[name]["argv"], cwd=root, stdout=fh, stderr=subprocess.STDOUT)  # noqa: S603
    _running[name] = proc
    if name.startswith(SHARED_KEY_PREFIX):
        take_key(name, f"registry job, pid {proc.pid}", pid=proc.pid)
    _meta_path(name).write_text(
        json.dumps({"pid": proc.pid, "started": time.time(), "finished": None, "code": None})
    )
    return status(name, root)


def status(name: str, root: Path = Path(".")) -> JobStatus:
    spec = CATALOG[name]
    meta = json.loads(_meta_path(name).read_text()) if _meta_path(name).exists() else {}
    proc = _running.get(name)
    if proc is not None and meta.get("pid") not in (None, proc.pid):
        # the job was started again elsewhere (CLI, another server): this process's handle is stale and
        # must not report that newer run as finished when the old one exits
        _running.pop(name, None)
        proc = None
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
    total_fn = spec.get("total_fn")
    if total_fn is not None:  # a total only knowable once the run has written its result
        spec = {**spec, "total": total_fn(root) or spec.get("total")}
    activity = last_activity(name, root)
    if state == "running" and activity is not None and activity > STALL_AFTER:
        state = "stalled"
    if state in ("unknown", "failed", "stalled") and is_complete(name, root):
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
        total = spec.get("total")
        if (
            state in ("idle", "unknown")
            and rp.exists()
            and (total is None or done < total)
            and time.time() - rp.stat().st_mtime < 600
        ):
            state = "running"
        if state in ("idle", "unknown") and total is not None and done >= total:
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
