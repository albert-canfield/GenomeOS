# Working on GenomeOS

## Branches

- `main` holds tested, working code only. It is protected: nothing lands there
  unless CI (lint + tests) is green.
- `dev` is where all work happens. Commit freely; CI runs on every push so
  breakage is visible early.
- When a piece of work on `dev` is complete and tested, open a pull request
  `dev -> main` (or run `scripts/promote.sh`), let CI pass, and merge.

```
git switch dev                 # work here
...edit, test...
uv run ruff check . && uv run ruff format --check . && uv run pytest -q
git commit -am "..."
git push
scripts/promote.sh             # PR dev -> main when green, then merge
```

## Every step

1. Tests for the change (`tests/`), real-data tests skip when the data is absent
   and use distilled summaries when possible (docs/STORAGE.md).
2. `uv run ruff check .` and `uv run ruff format .` clean.
3. Every rule or parameter added carries evidence and a confidence.
4. Record the outcome in docs/PROGRESS.md and, if something was learned from data,
   in docs/LESSONS.md.
