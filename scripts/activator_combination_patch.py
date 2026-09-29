# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Run anything under another activator-combination rule, and record where the rule matters.

Imported (as `activator_combination_patch`, with scripts/ on PYTHONPATH) it sets
`genomeos.runtime.grn.ACTIVATOR_COMBINATION` from GENOMEOS_COMB_RULE (default "mean"), makes
`runtime/located.py`, which repeats the mean formula inline, use the same rule, and records every
evaluation of a gene with two or more activators and whether the four candidate rules disagree
there. With GENOMEOS_COMB_SENS set, the record is written as JSON at exit, keyed by pytest node id.

    PYTHONPATH=scripts GENOMEOS_COMB_RULE=or GENOMEOS_COMB_SENS=/tmp/s.json \\
        uv run pytest -q -p activator_combination_patch tests
    PYTHONPATH=scripts GENOMEOS_COMB_RULE=max uv run python -c "import activator_combination_patch, \\
        sys; from genomeos.bio import main; sys.exit(main(['test', 'data/demo']))"

Scratch instrument for the census in scripts/activator_combination.py; it changes no committed value.
"""

from __future__ import annotations

import atexit
import inspect
import json
import os
import textwrap

import pytest

from genomeos.runtime import grn, located

RULE = os.environ.get("GENOMEOS_COMB_RULE", "mean")
if RULE not in grn.ACTIVATOR_RULES:
    raise ValueError(f"GENOMEOS_COMB_RULE={RULE!r}; one of {grn.ACTIVATOR_RULES}")
grn.ACTIVATOR_COMBINATION = RULE
_combine = grn.combine_activators
CURRENT = ["<outside a test>"]
SENS: dict[str, dict] = {}


def _record(terms: list[float], where: str) -> None:
    if len(terms) < 2:
        return
    vals = [_combine(terms, r) for r in grn.ACTIVATOR_RULES]
    d = SENS.setdefault(CURRENT[0], {"calls": 0, "rules_disagree": 0, "max_spread": 0.0, "where": set()})
    d["calls"] += 1
    d["where"].add(where)
    spread = max(vals) - min(vals)
    if spread > 1e-9:
        d["rules_disagree"] += 1
        d["max_spread"] = max(d["max_spread"], spread)


def _recording_combine(terms, rule=None):
    _record(terms, "grn")
    return _combine(terms, rule)


grn.combine_activators = _recording_combine

# located.py computes the mean inline; swap that expression for the shared function, verbatim otherwise
_INLINE = """            a = sum(
                r.strength * _hill(self._read(state, r.source, c), k[(r.id, c)], r.hill) for r in acts
            )
            a /= len(acts)"""
_SHARED = """            _t = [
                r.strength * _hill(self._read(state, r.source, c), k[(r.id, c)], r.hill) for r in acts
            ]
            _record(_t, "located")
            a = _combine(_t)"""
_src = textwrap.dedent(inspect.getsource(located.LocatedRuntime.fluxes))
if _INLINE not in _src:
    raise RuntimeError("runtime/located.py no longer computes the activator mean inline; update this patch")
_ns = dict(located.__dict__, _record=_record, _combine=_combine)
exec(_src.replace(_INLINE, _SHARED), _ns)  # noqa: S102 - the runtime's own source, one expression swapped
located.LocatedRuntime.fluxes = _ns["fluxes"]


def _dump() -> None:
    out = os.environ.get("GENOMEOS_COMB_SENS")
    if out:
        data = {k: {**v, "where": sorted(v["where"])} for k, v in SENS.items()}
        with open(out, "w") as fh:
            json.dump({"rule": RULE, "sensitive": data}, fh, indent=1, sort_keys=True)


atexit.register(_dump)


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_protocol(item, nextitem):  # noqa: ARG001 - pytest's signature
    CURRENT[0] = item.nodeid
    yield
    CURRENT[0] = "<outside a test>"
