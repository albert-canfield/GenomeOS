"""The panel the genome-wide roll-up normalises over may not shrink without saying so.

`normalise_family` turns each raw reading into a standardised residual of a fit over the panel of
biosamples, so the panel is that reading's denominator: the intercept, the slope and the standard
error are all computed over it, and dropping one biosample makes every residual in the file a
different number. `genome/reader.py` already refuses a panel with a missing *field* on the stated
ground that "a panel with a hole is a different panel". A missing *biosample* is the same hole, and
until `panel_or_refuse` nothing refused it, while two routes reached it in silence: `--rerun`
starts the summary afresh from the biosamples named on the command line and `DEFAULT_CELLS` is
narrower than the panel on disk, and the completeness filter drops any biosample absent from one
chromosome's row.
"""

import ast
import importlib.util
import json
import os

import pytest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_SRC = os.path.join(_ROOT, "scripts", "reader_genome_wide.py")
_SPEC = importlib.util.spec_from_file_location("reader_genome_wide", _SRC)
rgw = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(rgw)

PUBLISHED = [
    "K562",
    "HepG2",
    "GM12878",
    "H1",
    "IMR-90",
    "SK-N-SH",
    "cardiac muscle cell",
    "keratinocyte",
    "hepatocyte",
    "astrocyte",
    "CD14-positive monocyte",
    "testis",
    "ovary",
]


def test_a_rerun_on_the_default_cells_REFUSES_and_names_both_biosamples_it_would_drop():
    """Plant for the first guard: the exact command the roll-up offers today."""
    with pytest.raises(ValueError) as e:
        rgw.panel_or_refuse(list(rgw.DEFAULT_CELLS), PUBLISHED, "--rerun")
    msg = str(e.value)
    assert "testis" in msg and "ovary" in msg, msg
    # the refusal states both sizes, so the reader does not have to count the names
    assert "11 biosamples" in msg and "drop 2 of the 13" in msg, msg
    # and why it refuses rather than what it refuses
    assert "standardised by a fit over the panel" in msg, msg


def test_the_completeness_filter_dropping_one_biosample_REFUSES():
    """Plant for the second guard: no --rerun, one biosample missing from one chromosome's row."""
    complete = [c for c in PUBLISHED if c != "ovary"]
    with pytest.raises(ValueError) as e:
        rgw.panel_or_refuse(complete, PUBLISHED, "the roll-up")
    assert "ovary" in str(e.value) and "drop 1 of the 13" in str(e.value)


def test_the_guard_passes_the_panels_that_are_not_a_shrink():
    # the same panel, in any order
    rgw.panel_or_refuse(list(reversed(PUBLISHED)), PUBLISHED, "x")
    # a wider panel: a biosample added is not a biosample dropped
    rgw.panel_or_refuse([*PUBLISHED, "liver"], PUBLISHED, "x")
    # the first run, with no published panel to shrink
    rgw.panel_or_refuse(list(rgw.DEFAULT_CELLS), [], "x")


def test_the_footgun_the_guard_exists_for_is_live_in_this_checkout():
    """DEFAULT_CELLS really is narrower than the panel on disk, so the refusal is not hypothetical.

    `data/results/reader_genome_wide.json` is a tracked repository file: if it is absent that is a
    broken checkout and this raises, rather than skipping and reporting a pass it did not earn.
    """
    p = os.path.join(_ROOT, "data", "results", "reader_genome_wide.json")
    with open(p) as fh:  # a missing tracked file is an error here, never a skip
        on_disk = json.load(fh)["cell_types"]
    assert len(on_disk) == 13
    assert set(rgw.DEFAULT_CELLS) < set(on_disk)
    assert sorted(set(on_disk) - set(rgw.DEFAULT_CELLS)) == ["ovary", "testis"]


def _calls(fn: ast.FunctionDef, name: str) -> list[int]:
    # ast.walk yields no particular order, so the line numbers are sorted before they are compared
    return sorted(
        n.lineno
        for n in ast.walk(fn)
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id == name
    )


def test_both_guards_are_reached_and_the_rerun_one_fires_before_the_first_write():
    """A guard the roll-up never calls would pass every test above and protect nothing.

    The `--rerun` call has to come before any `save_result`, because the roll-up overwrites the
    summary as each chromosome finishes: a check after the loop would fire on a file already
    replaced.
    """
    with open(_SRC) as fh:
        src = fh.read()
    main = next(n for n in ast.parse(src).body if isinstance(n, ast.FunctionDef) and n.name == "main")
    guards = _calls(main, "panel_or_refuse")
    saves = _calls(main, "save_result")
    assert len(guards) == 2, guards
    assert saves, "the roll-up writes nothing, so there is nothing to guard"
    assert guards[0] < min(saves), (guards, saves)
    assert guards[1] > min(saves), (guards, saves)
