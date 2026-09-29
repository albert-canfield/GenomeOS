# SPDX-License-Identifier: AGPL-3.0-or-later
"""Application-side import sources for BioLang programs.

`register()` tells the engine's parser how to resolve `import protein:SYMBOL`: the block comes from the
packaged human proteome (genomeos.lib.proteome), so a program can pull a cited protein definition in
without the engine knowing where it lives.

The same function is declared in pyproject.toml as the `protein` entry point of the engine's
`biolang_import_resolvers` group, so `bio` finds it wherever GenomeOS is installed without `register()`
being called; a packaged engine with no GenomeOS beside it finds nothing and says so.
"""

from __future__ import annotations

from genomeos.lang.parser import IMPORT_RESOLVERS


def resolve_protein(symbol: str) -> str:
    from genomeos.lib.proteome import block

    text = block(symbol)
    if "not in the packaged proteome" in text:
        raise ValueError(f"protein {symbol!r} is not in the packaged proteome")
    return f"module protein.{symbol.upper()}\n{text}"


def register() -> None:
    IMPORT_RESOLVERS.setdefault("protein", resolve_protein)
