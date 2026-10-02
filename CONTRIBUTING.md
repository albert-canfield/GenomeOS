# Working on GenomeOS

Open source, MIT licence. These rules apply to people and to automated
sessions alike; the working rules for several sessions sharing one checkout
are in [docs/ROADMAP.md](docs/ROADMAP.md) §8 and the reasons in
[docs/LESSONS.md](docs/LESSONS.md).

## Licence of what you contribute

Contributions are accepted under the licence of the file you are changing
(inbound equals outbound): Apache 2.0 for the BioLang engine, AGPL-3.0-or-later
for the GenomeOS application ([LICENSING.md](LICENSING.md)). By contributing you
also agree that the copyright holder may offer your contribution under a
separate commercial licence; that is what keeps dual licensing possible. New
files carry the `SPDX-License-Identifier` of the part they belong to.

## Branches

- `main` holds tested, working code only. It moves to a `dev` sha whose
  latest CI `test` run is green, through `scripts/promote_main.sh` (below,
  "Release checks"). Its branch rule names the same check as required, but
  the rule does not bind admins until the owner says so.
  *(2026-09-29: the rule binds admins now; the owner turned on
  `enforce_admins`. With it on, a direct push of a sha whose only green
  `test` run was manually triggered was refused, so `main` moves through a
  pull request from a branch frozen at the chosen sha; see "Release checks".)*
- `dev` is where all work happens. A push to `dev` runs the pre-push hook,
  not CI; CI runs on the tip of `dev` once a day.
- **Pull requests to `main` are opened by Albert, by hand, when he decides a
  state of `dev` is worth promoting.** Nobody else, and no script or session,
  opens, updates or merges a pull request. `scripts/promote.sh` exists for
  Albert's manual use only. While a pull request is open, every push to `dev`
  runs CI once, so keep them short-lived.
  *(2026-09-29: `scripts/promote.sh` is retired; it refuses and points to
  `scripts/promote_main.sh`.)*

## Release checks (2026-09-28)

Item S7 of the second external review. What enforces what, before this date:

| Where | What ran | What it bound |
| --- | --- | --- |
| push to `dev` | pre-push hook: `scripts/check.sh` on the pushed sha in a clean worktree | only a clone that installed it; `GENOMEOS_SKIP_CHECK=1` or `--no-verify` skips it; it runs with this machine's data caches and `.venv` (`UV_NO_SYNC=1`), not a fresh install |
| pull request into `main` | CI on `opened`, `reopened`, `ready_for_review` | a later push to the pull request was never tested |
| daily, 06:00 UTC | CI on the tip of `dev` | nothing: red every day from 2026-09-22 to 2026-09-28 and no promotion waited on it |
| promotion to `main` | the pre-push hook (it checks pushes to `main` too); branch rule requiring `test` | the rule has `enforce_admins` off and sessions push as the owner's admin account |

Three findings behind the table. The scheduled run was testing `dev`, not
`main`: the repository's default branch is `dev` (repository settings, read
2026-09-28), and every scheduled run since 2026-09-16 recorded the sha that
was the tip of `dev` at the time, including 2026-09-22 (`cfa6d98`, when
`main` was `8767d2d`) and 2026-09-28 (`db04d1f`, when `main` was
`45a5aa6`). This checkout's `origin/HEAD` still points at `origin/main`,
because a clone sets it once; `git remote show origin` reports `dev`. The daily run
has been red for seven days on one test,
`tests/test_cancer_alterations.py::test_a_five_prime_partner_keeps_its_own_ectodomain`,
which reads caches git ignores (`data/knowledge/vep/`,
`data/knowledge/proteins/`) and lacks the `needs_caches` skip its neighbours
carry. In a clean worktree without those caches it fails here the same way
("no compiled protein definition for ALK") while its neighbour skips; the
pre-push hook links this machine's caches into its worktree, so it let the
same shas through. Neither of the last two moves of `main` had a green `test`
on its sha when it landed: `36f36d6` (2026-09-22) had no run yet and its
five later runs were red, and `45a5aa6` (2026-09-27) has never had one.

After this date:

- **Pull requests**: `synchronize` is a trigger, so every push to an open
  pull request into `main` runs CI; `concurrency` cancels the run a newer
  push supersedes.
- **Daily**: still the tip of `dev`, now asserted. The first step of `test`
  fails a scheduled run that is not on `refs/heads/dev`, so a change of
  default branch is reported instead of silently testing another branch.
  The checkout stays at `github.sha` rather than `ref: dev`: a run records
  `github.sha` as the sha it tested, and the gate and the branch rule both
  read that record. Checking out `dev` by name could test a newer sha than
  the one recorded (a push between trigger and checkout), or, if the default
  branch changed, test `dev` under `main`'s sha. `main` is not checked daily:
  it moves only by promotion, so a daily run would retest an unchanged sha.
- **Push to `main`**: CI runs on the promoted sha where it lands, once per
  promotion. It is a record, not a gate: it runs after `main` has moved, whatever
  route the sha took.
- **Promotion** is the gate. `main` moves only to a sha that is on
  `origin/dev`, contains `origin/main` (a fast-forward), and whose most
  recent completed `test` check run by GitHub Actions concluded `success`
  (not `skipped` or `neutral`, which the branch rule would accept).
  `scripts/promote_main.sh` checks exactly that and is a dry run unless
  given `--push`:

  ```
  scripts/promote_main.sh --latest-green     # the newest green sha between main and dev, dry run
  gh workflow run ci.yml --ref dev           # none green: test dev's tip now, then look again
  scripts/promote_main.sh --push SHA         # the coordinator, on the owner's go, once a day at most
  # 2026-09-29: --push is retired and refuses; the line above stays until its clean-up. Instead:
  scripts/promote_main.sh --prepare SHA      # on the owner's go: pushes branch promote-<first 7 of SHA> at SHA, nothing else
  scripts/promote_main.sh --verify SHA       # after the owner merges: main has SHA's file tree and SHA as an ancestor
  ```

  The push still goes through the pre-push hook, so the promoted sha is
  checked locally as well; the hook stays, as one check of two.

  *2026-09-29, beside the above: passing this gate does not establish that
  GitHub will accept the push. With `enforce_admins` on, `--push 9a59faf`
  passed every check here and GitHub refused it ("Required status check
  \"test\" is expected"), although a successful `test` check run by GitHub
  Actions existed on that exact sha, from a manually triggered run on `dev`.
  A successful manually triggered check is insufficient to establish GitHub
  eligibility. Until the gate reads what GitHub accepts, a promotion is a
  pull request from a branch frozen at the chosen sha (e.g.
  `promote-9a59faf`), whose checks run on the pull request and which the
  owner merges; the coordinator then verifies `main`'s file tree and
  ancestry against the chosen sha.*

  *2026-09-29, later, beside both: the gate now reads the event of the
  `test` run it relies on and reports two different things. "CI pre-check
  passed": that run concluded `success`. "Eligible for protected
  promotion": that run was also triggered by `pull_request` or `push`. A
  green `workflow_dispatch` or `schedule` run on `dev` is a pre-check only;
  it does not make the sha eligible, and eligibility comes from the pull
  request's own run. GitHub's documentation ("Troubleshooting required
  status checks", section "Checks from some workflow jobs are not
  evaluated") evaluates checks for pull requests and rulesets only from
  runs triggered by `push`, `pull_request`, `pull_request_review`,
  `pull_request_target`, `deployment` or `deployment_status`; that the same
  rule decided the refused direct push of `9a59faf` is an inference, not
  established. `--push` is retired and refuses. `--prepare SHA` makes the
  same checks, then pushes one branch, `promote-<first 7 of SHA>`, at
  exactly that sha (refusing if the branch exists at another sha), and
  prints the compare URL for the owner to open and merge the pull request,
  and the `--verify` command for after the merge. It never opens a pull
  request, never merges and never touches `main`. `--verify SHA` checks
  that `origin/main` has SHA's file tree and has SHA as an ancestor. The
  tests in `tests/test_promote_main.py` check the script's decisions and
  git operations against a local remote; GitHub's enforcement is not
  simulated.*

  *2026-09-29, later still, beside the note above (the reviewer's
  correction): a `pull_request` or `push` run is a qualifying CI event;
  protected promotion eligibility not established. The exact revision, the
  required checks and the protection rules decide, and GitHub's decision
  is authoritative; the same holds for the pull request's own run. The
  gate's output still says "eligible for protected promotion" on that
  branch, followed by a line saying it is a qualifying CI event only; the
  wording is held by the checkout guard until 2026-10-01 22:39 BST and is
  reworded after that. Direct pushes to `main` stay refused, and promotion
  stays the owner's: the owner opens and merges the pull request.*
- **Recommended branch rule for `main`** (the owner's setting; no session
  changes it): keep `test` from GitHub Actions as the required check and
  force pushes and deletions blocked, and turn on "Do not allow bypassing
  the above settings". GitHub then refuses, for everyone, a push of a sha
  without a passing `test`, and accepts a direct push of one that passed on
  `dev` ("After all required status checks pass, any commits must either be
  pushed to another branch and then merged or pushed directly to the
  protected branch", GitHub, *About protected branches*). Do not turn on
  "Require a pull request before merging": it forbids the fast-forward. The
  cost: while `dev` is red, `main` cannot move; `origin/dev` already keeps
  the day's work on GitHub.
- **When each change applies.** GitHub runs a workflow as it stands in the
  event's commit. Schedules read the default branch, which is `dev`, so the
  daily assertion applies from the first scheduled run after this reaches
  `origin/dev`; `synchronize` applies to a pull request whose head contains
  it; the push-to-`main` run first fires on the promotion that carries this
  commit to `main`.

## The cycle: finish, check, commit, push

1. **Finish the piece of work.** A feature with its tests, or a data job with
   its result file, or a document. Not a half-state.
2. **Check locally, the way CI does**: `scripts/check.sh` runs lint, format,
   the whole test suite and `bio test` on every program. It must print
   `check: green`. Lint and format only your own files when others are
   editing the checkout (`scripts/check.sh path/to/your_file.py ...`); the
   tests always run in full.
3. **Commit to `dev`** with a plain message that says what changed and what
   it proved. No attribution trailers. In a shared checkout commit through a
   private index with explicit paths (below), never `git add -A`.
4. **Push once per finished piece**, not per file. The `pre-push` hook in
   this checkout checks out the commit being pushed into a temporary
   worktree and runs `scripts/check.sh` there, so another session's
   half-edited files neither block nor pass your push; a red result refuses
   the push, and a red push is what produces the failure notifications. Set
   `GENOMEOS_SKIP_CHECK=1` only for a documentation-only push. The hook is
   local to the checkout (`.git/hooks/pre-push`); `scripts/install-hooks.sh`
   installs it in a fresh clone.
5. **Record**: a line in docs/PROGRESS.md (append-only) with the evidence,
   a lesson in docs/LESSONS.md if data taught one, a row in
   docs/DECISIONS.md if a decision was taken, the ROADMAP counts when a
   milestone moves.

```
scripts/check.sh                               # green before every commit
export GIT_INDEX_FILE=/tmp/idx-$$; git read-tree HEAD
git add path/to/your/files                     # explicit paths only
TREE=$(git write-tree); OLD=$(git rev-parse HEAD)
C=$(git commit-tree "$TREE" -p "$OLD" -m "What changed and what it proved")
git update-ref refs/heads/dev "$C" "$OLD"; unset GIT_INDEX_FILE
git read-tree HEAD                             # refresh the SHARED index, see below
git push origin dev                            # the hook checks again
```

**Refresh the shared index after every commit.** `update-ref` moves the branch
without touching `.git/index`, so a commit made through a private index leaves
the shared one behind by exactly that commit's contents. After a hundred such
commits the shared index is a hundred commits stale, `git status` in an editor
reports hundreds of phantom changes, and a plain `git commit` from it would
delete every file added since the drift began. One `git read-tree HEAD` after
each commit keeps it honest, costs nothing, and it is the reason this project
has had to repair that index three times in one day.

For a shared file (`cli.py`, `server.py`, `index.html`, `README.md`,
`PROGRESS.md`) build the staged copy from `git show HEAD:FILE` plus your own
hunk, `git hash-object -w` it and `git update-index --cacheinfo`; announce
the hunk to the other sessions first.

For a shared **markdown** file this is done for you, by heading:

```
uv run python scripts/stage_section.py docs/ATTRIBUTION.md "## My section" ["## the next heading"]
```

It takes the file as your private index already holds it, replaces only the
section whose heading you name with the working copy's, and writes the result
into that index. A section runs to the next heading of the same level or
shallower, so `"### E. From one cell to an organism"` takes area E and not the
rest of the roadmap; headings inside fenced code blocks are ignored. Nothing outside your section can reach the commit, because
everything outside it is the base verbatim; run it once per section you own and
the runs compose. It refuses, staging nothing, when your heading is missing from
the working copy or when a new section has no named heading to go before, rather
than appending to the end of a file another lane is also writing. This exists
because `git add docs/ROADMAP.md` has twice published a lane's half-written text
under another lane's commit.

## Every change

1. Tests for the change (`tests/`); real-data tests skip when the data is
   absent and read the distilled summaries when they exist (docs/STORAGE.md).
2. Every rule or parameter added carries evidence and a confidence;
   `UNKNOWN` stays `UNKNOWN`.
3. Raw downloads never enter the repository: stream, distil, discard.
4. Long jobs run once, through the jobs registry, and commit one result
   file per chromosome as it lands.
5. The engine (`genomeos/lang`, `ir`, `runtime`, `std`, `bio.py`) imports
   nothing from the genome, knowledge or web layers.

## A result written in this shared checkout (2026-10-02)

Several sessions work in one tree, so a lane's result manifest picks up every
other lane's uncommitted files and `result_manifest.code.dirty` is true through
no fault of its own. Requiring `dirty: false` outright would force the lanes to
run one at a time, so the rule is this instead. A result whose stamp is dirty is
sound only when it also carries, computed rather than asserted:

1. its **own** code committed (`own_uncommitted_code` empty), with the sha;
2. the **foreign** uncommitted files named;
3. the **transitive import closure** of its entry script, **computed by code**,
   showing none of those foreign files on it. A hand-written list of imports
   does not meet this: the closure is what decides whether a dirty file could
   have entered the number, and only code can enumerate it.

And before any number from such a result is quoted in README, docs/ROADMAP.md,
docs/ATTRIBUTION.md or a reply to the owner, run `scripts/manifest_rebuild.py`
on it in a clean worktree at its own sha.

**What the rebuild must show (amended 2026-10-02).** `0` differences **at every
leaf**, with one class of exception: the fields that describe the tree the run
happened in, by **exact path** and never by pattern — `dirty` and the foreign,
dirty and untracked path lists, under `code_cleanliness` and under
`result_manifest.code`. They are listed in the report under
`environment_fields_ignored`, and the exact paths are `ENVIRONMENT_FIELDS` in
the script rather than a sentence here, so the rule is a mechanism.

The first wording asked for 0 differences full stop, and that could only ever be
failed by a result honest about a shared checkout: it records the other lanes'
outstanding files, and a clean worktree has none, so the honest result failed
while one that omitted the record passed.

Two fields must hold on **both** sides, so the exception can never excuse a
result no commit reproduces: `own_code_is_committed` is true, and
`foreign_uncommitted_code_on_the_counting_path` is empty. Any other difference
is a defect.

Report **leaves compared, not top-level keys**. A key count hides a nested
difference: `data/results/cell2_eligibility.json` has 14 top-level keys and
**707** leaves. A result that fails any of this is a defect, not a result.
