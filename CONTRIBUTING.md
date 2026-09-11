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

- `main` holds tested, working code only. It is protected: nothing lands there
  unless CI (lint + tests + self-testing programs) is green.
- `dev` is where all work happens. Every push to `dev` runs CI.
- **Pull requests to `main` are opened by Albert, by hand, when he decides a
  state of `dev` is worth promoting.** Nobody else, and no script or session,
  opens, updates or merges a pull request. `scripts/promote.sh` exists for
  Albert's manual use only. An open pull request makes every push to `dev`
  run CI twice, so keep them short-lived.

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
git push origin dev                            # the hook checks again
```

For a shared file (`cli.py`, `server.py`, `index.html`, `README.md`,
`PROGRESS.md`) build the staged copy from `git show HEAD:FILE` plus your own
hunk, `git hash-object -w` it and `git update-index --cacheinfo`; announce
the hunk to the other sessions first.

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
