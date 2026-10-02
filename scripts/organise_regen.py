# SPDX-License-Identifier: AGPL-3.0-or-later
"""Regenerate one chromosome's `organised_<chrom>` result in place, through the writer that owns its
declaration (`genomeos.attribution.organise.run_and_save`, strict).

This script exists so that the result's own `code.argv` names something a second checkout can run.
The 24 committed `organised_chr*` results were written with no manifest at all, so they declared no
command; regenerating them from an ad-hoc script outside the repository would have replaced that gap
with a declaration that cannot be followed -- `scripts/manifest_rebuild.py` runs `uv run python` over
`code.argv` inside a worktree of the recorded `code.git_sha`, so an argv[0] outside the checkout is a
path that is simply not there. Twenty-five other results in the registry record `scripts/<name>.py`;
this is the same shape for this one.

`--top` IS REQUIRED and has no default, which is the point of having it on the command line at all.
The committed results are heterogeneous: chr21 was written with top=3 and the other 23 with top=25,
while the `organise` subcommand defaults to 15 and the function to 25, so NO committed file matches
the CLI default and a regeneration that accepted a default would silently resize chr21's lists (3
candidates to 5, 3 largest copies to 25) and truncate nothing on the rest only by luck. The value is
recorded in the result manifest's `parameters.top`, where the committed files recorded it nowhere.

    uv run python scripts/organise_regen.py --chrom chr21 --top 3
    uv run python scripts/organise_regen.py --chrom chr1 --top 25

Nothing is read under data/ before the writer runs, so the input-tracing window `genomeos.manifest`
arms at import holds exactly the reads the writer made (`result_manifest.traced_inputs`).
"""

from __future__ import annotations

import argparse
import sys


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--chrom", required=True, help="one chromosome, e.g. chr21")
    ap.add_argument(
        "--top",
        required=True,
        type=int,
        help="candidates and largest copies kept in the result file; no default on purpose",
    )
    args = ap.parse_args(argv)
    if args.top < 1:
        ap.error("--top must be at least 1")

    from genomeos.attribution.organise import run_and_save

    out = run_and_save(args.chrom, top=args.top)
    ru = out["real_unknown"]
    print(
        f"{args.chrom}: {out['blocks']} UNKNOWN blocks, {out['copies']} copies, "
        f"{ru['blocks']} real unknown ({ru['bp'] / 1e3:.0f} kb), "
        f"{len(out['candidates'])} candidates and {len(out['largest_copies'])} largest copies at "
        f"top={args.top}; {len(out['inputs']['declared'])} inputs declared, "
        f"budget branch {out['inputs']['budget_branch']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
