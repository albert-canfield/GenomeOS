# SPDX-License-Identifier: AGPL-3.0-or-later
"""The coordination screen: the two collisions worth reporting, and the staleness line."""

from __future__ import annotations

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("coordinate", ROOT / "scripts" / "coordinate.py")
assert spec and spec.loader
co = importlib.util.module_from_spec(spec)
spec.loader.exec_module(co)


def entry(who: str, hours: float, files: list[str], state: str = "running") -> dict:
    return {"who": who, "age": hours * 3600.0, "files": files, "state": state, "task": "t", "area": "J"}


def test_a_file_two_sessions_declare_is_reported_as_a_collision() -> None:
    """Two lanes editing one file is the merge this project has already been bitten by."""
    board = [
        entry("genomeos-aa", 0.5, ["genomeos/compare.py"]),
        entry("genomeos-bb", 0.5, ["genomeos/compare.py"]),
    ]

    out = co.collisions(board, [])

    assert out["files_claimed_by_more_than_one_session"] == {
        "genomeos/compare.py": ["genomeos-aa", "genomeos-bb"]
    }


def test_a_dirty_file_no_live_session_declared_is_reported() -> None:
    """The four files left dirty overnight were exactly this: owned by sessions that had stopped.

    Committing one guesses on the owner's behalf, which is how a twin lost its provenance.
    """
    board = [entry("genomeos-old", 40.0, ["data/twins/HG002.json"])]  # stale, so not live

    out = co.collisions(board, [{"path": "data/twins/HG002.json"}])

    assert out["dirty_files_no_live_session_has_declared"] == ["data/twins/HG002.json"]


def test_a_dirty_file_its_live_owner_declared_is_not_a_collision() -> None:
    """A lane editing what it said it would edit is the normal case and must stay quiet."""
    board = [entry("genomeos-now", 0.2, ["genomeos/evidence.py"])]

    out = co.collisions(board, [{"path": "genomeos/evidence.py"}])

    assert out["dirty_files_no_live_session_has_declared"] == []
    assert out["files_claimed_by_more_than_one_session"] == {}


def test_one_session_declaring_a_file_twice_is_not_a_collision() -> None:
    """Its own earlier entry is not a peer; only distinct sessions collide."""
    board = [entry("genomeos-aa", 0.5, ["a.py"]), entry("genomeos-aa", 30.0, ["a.py"])]

    assert co.collisions(board, [])["files_claimed_by_more_than_one_session"] == {}
