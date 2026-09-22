# SPDX-License-Identifier: AGPL-3.0-or-later
"""Application-side import sources for BioLang programs.

`register()` tells the engine's parser how to resolve `import protein:SYMBOL`: the block comes from the
packaged human proteome (genomeos.lib.proteome), so a program can pull a cited protein definition in
without the engine knowing where it lives.
"""

from __future__ import annotations

from genomeos.lang.parser import IMPORT_RESOLVERS


def _protein(symbol: str) -> str:
    from genomeos.lib.proteome import block

    text = block(symbol)
    if "not in the packaged proteome" in text:
        raise ValueError(f"protein {symbol!r} is not in the packaged proteome")
    return f"module protein.{symbol.upper()}\n{text}"


def register() -> None:
    IMPORT_RESOLVERS.setdefault("protein", _protein)
