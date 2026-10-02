# False red verdicts, 2026-10-02

A red verdict in `.git/genomeos-check/` is a statement about a tree. Four of the reds written on
2026-10-02 were statements about nothing in any tree: two because `git check-ignore` could not answer
past a symlinked data store, and two because `ruff` was handed a `.json` as a lint argument. They were
read as real, one of them refused a push, and `origin/dev` sat thirty commits behind for an hour.

This file exists because a later reader cannot tell a tooling-caused red from a test red by looking at
one, and because the status files themselves do not survive: they are pruned after a fortnight and are
keyed by tree hash alone, so a peer who runs the check on the same tree overwrites one. Both of those
things happened while this file was being written, and both are recorded below rather than smoothed
over.

The mechanisms that stop the recurrence are elsewhere and are not notes:

- `tests/local_data.py` answers a symlinked store BY NAME instead of raising
  (`tests/test_local_data_symlinked_store.py`).
- `scripts/check.sh` refuses a non-`.py` lint argument by name and does not kill the run
  (`tests/test_check_sh_arguments.py`).
- every status file now carries `error_class`, so the distinction below is in the record rather than
  in this file (`genomeos/verdict.py`).
- `.git/hooks/pre-push` is a fixed wrapper that runs the COMMITTED checker, so the stale copy that
  caused the first two cannot recur (`scripts/pre-push-hook.sh`,
  `tests/test_pre_push_hook_is_the_wrapper.py`).

## The two symlink reds

Both are in the record, both say 20 errors and 0 failures, and neither is about any code.

| tree | finished | counts | where it ran |
| --- | --- | --- | --- |
| `7446015a503fcc6f281173590591709d1f6e7aab` | 18:13:43 | 4,418 passed, 0 failed, **20 errors**, 19 skipped | a lane's worktree, `scratchpad/wt6` |
| `579b4949ea024ef98722a9ea14f9825478818724` | 18:53:59 | 4,436 passed, 0 failed, **20 errors**, 11 skipped | the push's worktree, `genomeos-push.A2u4CX` |

Every one of the forty errors was the same line, on `data/cache/entex/alphagenome_track_metadata_copy.csv`:

```
local_data.NotMachineLocalError: git could not say whether <path> is ignored
(fatal: pathspec '<path>' is beyond a symbolic link)
```

`git check-ignore -q --no-index <path>` exits **128** when the path is reached through a symlink and
**0** against a real directory. `tests/local_data.py` treated any code but 0 or 1 as unanswerable and
raised, correctly by its own rule, so twenty tests ERRORED at setup instead of running.

**The cause was a stale hook, and two other diagnoses were wrong.** `.git/hooks/pre-push` was a *copy*
taken at 04:34 which still did `ln -s "$root/data/$d" "$tmp/data/$d"`, while `scripts/pre-push.sh` had
moved at 17:52 to read-only APFS clones. `scripts/install-hooks.sh` copies, so the hook was not the
script and the clone switch never reached what the push ran. Measured refutations, kept because a fix
aimed at either would have passed its own test while every push stayed red:

- **`$TMPDIR` behind `/var -> private/var` is not the cause.** A worktree under `$TMPDIR` whose stores
  are real directories answers 0 — measured with the main repository also under `/var` and with it on
  an unsymlinked path. git resolves a linked worktree's root physically and a relative pathspec is
  taken from the process's physical cwd.
- **`realpath` before `check-ignore` does not fix it**, and fails in the dangerous direction. An
  ignore rule is about a *name*, so resolving the path asks a different question. In the exact failing
  shape the resolved path lies in another working tree of the same repository and git exits 128 again
  (`is outside repository at '<worktree>'`); where the target sits inside the same worktree under
  another name, git exits **1**, "not ignored", and a correct marker is reported as wrong. `--no-index`
  changes neither reading.

Measured before and after, in a worktree of `d55534c` with the three stores symlinked exactly as the
stale hook did them: the five files that held the twenty errors go from **177 passed / 20 errors** to
**197 passed / 0 errors**. The tests RUN, because the stores are present through the link. A store that
is genuinely absent behind a symlink skips with the store named and the fetch command given.

## The ruff-on-JSON reds: one in the record, one observed and no longer there

`ruff` is a Python linter and reads whatever path it is given as Python. A lane passed
`data/results/manifest_headlines.json`; ruff reported **96 errors** in it, `set -e` stopped the run at
the `ruff-check` leg, no test ran at all, and a red file was written for the tree.

| tree | finished | verdict | counts | scope files |
| --- | --- | --- | --- | --- |
| `b17ef4d0b2800f21f64e34c4dbee62ea707cc180` | 17:58:06 | red, leg `ruff-check` | **`null`** — no test ran | `genomeos/cancer/tumour.py`, `scripts/constrained_unknown_targets.py`, `data/results/constrained_unknown_targets_v2.json`, `data/results/manifest_headlines.json` |

`counts: null` is the signature: the run died before pytest, so the file carries no evidence about the
suite in either direction.

**The second one is not in the record, and that is worth more than the second row would have been.**
This lane was told there were two and read `status-f8a9cf884bdf2f665f9dde009400d718749dc560` as
`red`, leg `ruff-check`, scope naming `data/results/manifest_headlines.json` and
`tests/test_headline_registry_matches_readme.py`. On re-reading minutes later that path held a
different run — `tree_moved`, 4,497 passed, 0 failed, 0 errors, scope naming only the test file, mtime
18:58:16. A sweep of every status file in the directory now finds **exactly one** whose `scope_files`
contain a `.json`. The status file is keyed by tree hash alone, so a peer checking the same tree
overwrites it; the two readings cannot both be of the same bytes, and this file does not claim to know
which was which. What it records is that the count of this class in the record is **one**, that a
second was seen, and that a verdict store keyed by tree alone cannot be used to count anything over
time. Neither reading changes the mechanism or the fix.

## The red that was a peer saving a file, and the only one with a `tree_moved` mechanism

`tests/test_context_evidence.py::test_the_census_script_s_own_closure_is_the_module_s_closure`
failed once at 22:07 in a run of eleven test modules: **305 passed, 1 failed.** The assertion that
broke compared two `foreign_uncommitted_code_on_the_counting_path` lists and the diff named
`genomeos/predict/alphagenome_adapter.py` — a file **the test does not mention**, belonging to a
lane that was editing it at that moment.

The test called `code_cleanliness` **twice** and asserted the two results agree, once with the
census entry spelled relatively and once absolutely. Each call runs its own `git status`. In a
checkout several sessions share, a peer saving a file between the two makes them disagree, and the
verdict is then a statement about **no tree either call measured** — the `tree_moved` class, which
until now had no worked example in this file. It passed three times out of three on re-run, which
is the signature and is also why such a red is tempting to dismiss rather than explain.

**The mechanism was planted rather than inferred**, because "a peer must have saved something" is a
story and three green re-runs are not evidence of a cause. In a planted repository with one
committed script importing one committed module, two consecutive `code_cleanliness` calls with the
module modified between them return:

| | `counting_path` | `foreign_uncommitted_code_on_the_counting_path` |
| --- | --- | --- |
| first call | identical | `[]` |
| second call | identical | `['genomeos/peer.py']` |

The closure is **not** what moved — it comes from the imports — and that is exactly the shape the
one red had, with the first assertion passing and the second failing.

**The fix is that both readings now come from one snapshot:** `code_revision` is called once and
both calls are given its answer, so the only thing that differs between them is the entry's
spelling, which is what the test was ever about. The plant became the regression test
(`test_two_readings_of_one_tree_disagree_when_it_moves_between_them`), so the fix carries a reason
that can be checked instead of a comment asserting one.

**The general form, and it is this file's fourth instance of the same shape:** a check that reads
the live checkout twice and compares the readings is not testing the code, it is testing whether
anyone else was working. `tests/test_code_cleanliness_hermetic.py` exists because of the same
defect one step further out, and `4f44dbf` fixed a third where a historical premise was asserted
against the live tree. A claim about a tree belongs to that tree, by sha or by planting.

## Two reds the coordinator bought by choosing WHEN instead of WHAT (2026-10-02)

Both are pushes, both are the coordinator's, and both produced verdicts about trees nobody
intended to ship. They are the same mistake twice: **the sha was decided by when something ran
rather than by what was judged.**

**The mid-sequence push.** push17 ran on a sha taken from `HEAD` while a lane was between two of
its own commits. A lane that commits code and then regenerates the result that code governs passes
through a tree where the new code and the old result disagree — here the module's `FALSIFIER`
already carried amendment 2's wording while the committed registration still carried amendment 1's,
because the result landed two commits later. **Five tests failed, and none of them was about the
code.** A tree between a lane's code commit and its result commit is red BY CONSTRUCTION.

**The late sha capture.** push18 was deliberately queued behind a timed wait, so the
one-push-per-thirty-minutes interval would be measured from the previous push's start. It was
launched as `--sha $(git rev-parse HEAD)`. **That substitution evaluated when the WAIT fired, not
when the decision was made**, and a peer had committed twice in between — so the announced sha and
the pushed sha were different (`dbb5d4a` announced, `47835d1` pushed). It happened to land on a
boundary: the two extra commits were one lane's registration-only sequence, `dbb5d4a` was an
ancestor, and all 52 files of the closure being verified had identical blobs at both shas, checked
independently. **That is luck, not method.** A push is pinned to a sha precisely so the tree cannot
move under the verdict, and the movement was reintroduced inside the command that does the pinning.

**The rule, in two parts:**

- **Push a sha at a SEQUENCE BOUNDARY** — a lane's declared-green sha, or `HEAD` when no lane is
  mid-sequence — and ASK rather than assume. A registration-only commit is a boundary by
  construction, because there is no result for it to disagree with.
- **Capture the sha at DECISION time and pass it literally.** Never let a command substitution
  resolve the subject of a verdict at execution time.

**Why this belongs in this file rather than in a runbook:** neither red said anything about the
code, both looked exactly like real failures, and the first one cost an hour of classification work
before anybody noticed the tree was the problem. That is the definition this file exists for.

## Telling them apart from now on

`genomeos/verdict.py` writes an `error_class` into every status file, with the list of classes beside
it so a reader needs neither this file nor the source:

| class | meaning |
| --- | --- |
| `test` | a test FAILED. Something the suite asserts about the code did not hold. |
| `test_setup` | pytest red with ERRORS and NO failures: tests could not START. Where environment and tooling faults land — **both symlink reds above are this**. |
| `tooling` | a non-test leg refused or failed: lint, format, shellcheck, startup. **The ruff-on-JSON red is this.** |
| `tree_moved` | the working tree changed mid-run, so the counts belong to no one tree. |
| `unknown` | red or unreadable with nothing in the file to attribute it to. Never a pass. |

It is a reading aid and nothing more. `verdict require` still refuses anything that is not green,
whatever class it carries, because a red verdict of any class cannot answer for a push.

## What this file does not say

It does not say these were the only false reds of the day. The coordinator's reading is that any
verdict taken before about 19:00 may have been produced under the stale symlinking hook, and nobody
has gone back through them. It also does not say that a red is usually tooling: four reds in one day
were, and the rest were real.
