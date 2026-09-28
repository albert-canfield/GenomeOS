# SPDX-License-Identifier: Apache-2.0
"""A key a block does not read is an error that names it, not a silent no-op.

Until 2026-09-27 the parser stored every `key: value` and read only the ones it knew, so the engine's
smoke program ran for months with `basal_rate:` (the key is `basal`) and `rate:` on its rules (a rule's
weight is `strength`), and passed on an all-zero run.
"""

import pytest

from genomeos.lang import BioLangError, grammar, parse
from genomeos.lang import parser as lang_parser


def _error(src: str) -> str:
    with pytest.raises(BioLangError) as e:
        parse("module x\n" + src)
    return str(e.value)


def test_unknown_gene_key_names_key_block_accepted_keys_and_closest():
    msg = _error("gene A { basal_rate: 1.0 }")
    assert "gene 'A' has no key 'basal_rate'" in msg
    assert "did you mean 'basal'?" in msg
    assert "A gene accepts: " in msg and "basal" in msg and "evidence" in msg


def test_unknown_key_without_a_close_match_still_lists_the_accepted_keys():
    msg = _error("gene A { basal: 1 }\nprotein Ap { }\nrule A produces Ap { rate: 1.0 }")
    assert "rule 'A produces Ap' has no key 'rate'" in msg and "did you mean" not in msg
    assert "strength" in msg


def test_nested_transcript_keys_are_checked():
    assert "did you mean 'exons'?" in _error("gene A { transcript t { exon: chr1:1-2 } }")


def test_every_block_kind_accepts_exactly_the_grammar_table():
    for kind in lang_parser._KINDS:
        assert set(lang_parser._accepted(kind)) == set(grammar.BLOCKS[kind].get("props", {})) | set(
            grammar.COMMON
        )


def test_a_transcript_without_evidence_carries_its_genes():
    m = parse(
        'module x\ngene A { transcript t { exons: chr1:1-20 }\n evidence: curated "src"; confidence: 0.9 }'
    )
    t = m.entities["t"]
    assert t.evidence.kind.value == "curated" and t.evidence.source == "src" and t.confidence == 0.9


def test_a_transcript_with_its_own_evidence_keeps_it():
    m = parse(
        'module x\ngene A { transcript t { exons: chr1:1-20; evidence: predicted "m" }\n'
        ' evidence: curated "src"; confidence: 0.9 }'
    )
    assert m.entities["t"].evidence.kind.value == "predicted"
