# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Whether THIS ENVIRONMENT provides a top-level module, asked of the finders and not of `sys.modules`.

One function, because the project had the same question answered two different ways and one of the two
answers could be moved by anything that had put a module object in `sys.modules`. The gate that decides
whether an optional, paid client is available is not a place for a second opinion.

`importlib.util.find_spec`, which both call sites used until 2026-10-03, consults `sys.modules` BEFORE
it asks the filesystem, and it RAISES rather than answering when what it finds there has no spec:

    >>> sys.modules["alphagenome"] = types.ModuleType("alphagenome")   # __spec__ is None
    >>> importlib.util.find_spec("alphagenome")
    ValueError: alphagenome.__spec__ is None

A `types.ModuleType` stand-in is exactly what tests/test_adapter_gene_axis.py:136 and
tests/test_model_version_pin.py:36-50 build to stand in for the uninstalled client, and measured on
this machine with the real client INSTALLED, `predict.status()` raised that ValueError instead of
reporting `package: True`. A stub carrying a real spec moves the answer the other way: an absent module
reads as present, which is how a test that must skip runs against a stand-in instead.

The finders are asked in `sys.meta_path` order, which is what an actual `import` consults AFTER the
`sys.modules` shortcut, so a module provided by an editable install's own finder still answers True
while nothing in `sys.modules` can change the answer. `importlib.machinery.PathFinder` alone would be
narrower but would miss editable installs.

The residual blind spot, named rather than closed: code that installs a `MetaPathFinder` of its own and
leaves it on `sys.meta_path` moves this answer too. That is a deliberate import hook rather than a
leftover module object, and nothing in this project installs one.

This answers about the environment and NOT about whether the module imports cleanly: a module whose
import raises is still importable by this reading. The call sites want to know whether an optional
extra is installed, which is what a finder answers; whether it then works is the import's own business.
"""

from __future__ import annotations

import sys

#: What a finder may raise about a name it cannot resolve. An unaskable question is answered False and
#: not raised, because a gate that throws is not a gate that says no: `find_spec` raising `ValueError`
#: through `predict.status()` is the whole reason this module exists.
_FINDER_ERRORS = (ImportError, AttributeError, TypeError, ValueError)


def importable(module: str) -> bool:
    """True if some finder on `sys.meta_path` can locate `module` in this environment."""
    for finder in sys.meta_path:
        try:
            spec = finder.find_spec(module, None)
        except _FINDER_ERRORS:
            continue
        if spec is not None:
            return True
    return False
