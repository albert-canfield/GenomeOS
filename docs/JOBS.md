# Background jobs: reliable by construction

**Status on 2026-10-10.** Accurate as written, checked against `genomeos/jobs.py`: the pid
refusal, the heartbeat, the 20-minute stall threshold (`STALL_AFTER`), the supervisor's
auto-heal and the append-mode logs are all in the registry, and completion is still a
predicate over the results (`complete`), not over the process — the proteome job's is still
25 `proteome_chr*.json` files. The registry has grown since this page was written, so the
list of jobs lives in the code and in the Progress tab rather than here.

Long computations (a genome-wide pass, a proteome compile) run as separate
processes started from the Progress tab or the command line and tracked in
`data/jobs/`. The rules that make them reliable:

| rule | mechanism |
|---|---|
| resumable | every job saves per unit (a chromosome, a gene) into `data/results` or the knowledge cache and skips what is already saved on the next start |
| never gives up | the genome-wide scripts loop over passes with a growing pause (60 s to 30 min) until every unit has its result; a source outage costs time, not the job |
| one copy only | the registry records the pid on disk and refuses to start a twin while that pid is alive, even after the server restarts |
| proof of life | jobs write a heartbeat (`data/jobs/<name>.heartbeat`) at every step; the Jobs panel shows "active N s ago" and the position inside the current step (proteome: `chr19: 320/1,397`) |
| stall detection | a running job with no heartbeat or log line for 20 minutes is shown as `stalled` |
| self-healing | the server's supervisor thread checks every minute: an auto-heal job that is not complete and is dead or stalled is restarted from what it saved, at most once per 10 minutes; the panel counts heals |
| manual nudge | the **Heal** button on a failed or stalled job does the same on demand; `force` restarts a live job |
| clean logs | logs are opened in append mode, so a log can never show an older process's lines after the live ones |

Completion is a fact about the results, not the process: the proteome job
is complete when 25 chromosome coverage files exist, the UNKNOWN job when
24 chromosomes carry the curated repeat pass. A job that is complete is
never restarted.
