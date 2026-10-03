# SPDX-License-Identifier: AGPL-3.0-or-later
"""Version an already-published compiled program whose header sentence was false.

The paragraph quarantined here as `CELL_PARAGRAPH_SUPERSEDED` claimed that a rule's
`when: cell_type` names a cell the element had been assayed or scored within. For a predicted rule
it does not: the compiled cell is the argmax over the scorer's RNA-seq track axis, the tissue whose
predicted expression of the TARGET GENE moved most on deletion. On chr21 that is 5,174 of 5,176
rules, and exactly one shipped program carries the paragraph.

This module is the ONLY place in the package that holds a copy of that text, so a test can assert
the compiler itself is free of it and that no compiler run can emit it again.

The repair is VERSIONED, never in place. `data/results/label_gene.json` and
`data/results/label_gene_registration.json` name `data/organisms/human/noncoding_chr21.bio` by
sha256, so v1's bytes are a pin and rewriting them would break a rebuild.
"""

from __future__ import annotations

from genomeos.attribution.compile import cell_sentence_lines, count_rules

#: The superseded sentence, verbatim as the compiler wrote it and as v1 of the program carries it.
#: Quarantined: nothing builds program text from it, and a test asserts it is absent from the
#: compiler and from the new program.
CELL_SENTENCE_SUPERSEDED = (
    "Every rule is gated on the cell it was measured or predicted in (`when: cell_type = K562`), one"
)

#: The whole superseded paragraph as v1 of a compiled program carries it, `# ` included.
CELL_PARAGRAPH_SUPERSEDED = (
    "# " + CELL_SENTENCE_SUPERSEDED,
    "# rule per element, gene and cell, so a run in HepG2 integrates none of K562's. A rule whose",
    "# cell was not recorded says `cell_type = unknown`, which matches no cell: it is never universal.",
)


def reheader_cell_sentence(text: str) -> str:
    """Return `text` with the superseded cell paragraph replaced, and NOTHING else changed.

    The two counts in the new paragraph are read off `text`'s own rule lines, read off the result
    again and compared; and the result's rule lines must be byte-identical to the input's. This
    refuses on either, so it cannot write a header that disagrees with its own rules and cannot
    move a figure in a published program while repairing its prose.
    """
    lines = text.splitlines()
    para = list(CELL_PARAGRAPH_SUPERSEDED)
    at = next((i for i in range(len(lines) - len(para) + 1) if lines[i : i + len(para)] == para), None)
    if at is None:
        raise ValueError("the superseded cell paragraph is not in this program, verbatim")
    predicted, experimental = count_rules(text)
    out = lines[:at] + cell_sentence_lines(predicted, experimental) + lines[at + len(para) :]
    result = "\n".join(out) + ("\n" if text.endswith("\n") else "")
    if count_rules(result) != (predicted, experimental):
        raise ValueError("the rewritten program's rule counts differ from the ones written into it")
    before = [ln for ln in lines if ln.startswith("rule ")]
    after = [ln for ln in out if ln.startswith("rule ")]
    if before != after:
        moved = sum(1 for a, b in zip(before, after, strict=False) if a != b)
        raise ValueError(
            f"{moved} rule lines moved and {len(after) - len(before):+d} were added or dropped; "
            "a header fix changes no rule"
        )
    if any(CELL_SENTENCE_SUPERSEDED in ln for ln in out):
        raise ValueError("the superseded sentence is still in the rewritten program")
    return result
