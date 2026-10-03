#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Write a new version of a compiled BioLang program with its cell paragraph repaired.

The superseded paragraph said each rule's cell was one the element had been assayed or scored
within. For a predicted rule that is false: the compiled cell is the argmax over the scorer's RNA-seq track
axis, so the element was predicted in every track and in no one of them. On chr21 that is 5,174 of
5,176 rules.

The repair is VERSIONED, never in place. `data/results/label_gene.json` and
`data/results/label_gene_registration.json` name `data/organisms/human/noncoding_chr21.bio` by
sha256, so v1's bytes are a pin and rewriting them would break a rebuild. This reads v1, replaces
the paragraph, and writes a new path; v1 is opened read-only.

It also REFUSES rather than ship a header that disagrees with its own body: the two counts in the
new sentence are read off the input's rule lines, read off the output again and compared, and the
output's rule lines must be byte-identical to the input's. See
`genomeos.attribution.reheader.reheader_cell_sentence`.

    uv run python scripts/reheader_program.py <in.bio> <out.bio>
"""

from __future__ import annotations

import sys
from pathlib import Path

from genomeos.attribution.compile import count_rules
from genomeos.attribution.reheader import reheader_cell_sentence


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print(__doc__, file=sys.stderr)
        return 2
    src, dst = Path(argv[1]), Path(argv[2])
    if src.resolve() == dst.resolve():
        print(f"refusing to rewrite {src} in place: its bytes are a pin", file=sys.stderr)
        return 1
    text = src.read_text()
    out = reheader_cell_sentence(text)
    predicted, experimental = count_rules(out)
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(out)
    print(
        f"wrote {dst} ({len(out.encode())} bytes) from {src}: "
        f"{predicted} predicted and {experimental} experimental rule lines, both counted, "
        "rule lines byte-identical to the input"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
