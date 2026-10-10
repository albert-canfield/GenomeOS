# SPDX-License-Identifier: AGPL-3.0-or-later
"""Every shipped `design` program, run through BioForge, with its answer recorded so CI can hold it stable.

Area H, 2026-10-10. `bio test` runs a program's `experiment` blocks and asserts, not its `design` blocks,
so a change to a module could move a design's answer without any check noticing. This module finds every
design block in the paths CI's `bio test` step covers, runs each once at a fixed seed (cached, so the
tests that read the shipped designs share one run per process) and records the answer: the best
feasible perturbation, the candidates tied with it, the organism runs spent and the feasible candidates
found. tests/test_forge_design_answers.py compares the record with a pinned one.

The record states what BioForge answers, never whether the answer is right: no field in it is an
accuracy, a match to a publication or a probability (NOT_A_CLAIM).
"""

from __future__ import annotations

import re
from functools import cache
from pathlib import Path
from typing import Any

from genomeos.lang import parse_file
from genomeos.organism.forge import DesignResult, run_designs

# the paths `uv run bio test data/demo genomeos/std data/organisms` covers in .github/workflows/ci.yml
PROGRAM_ROOTS = ("data/demo", "genomeos/std", "data/organisms")
DESIGN_BLOCK = re.compile(r"^[ \t]*design[ \t]+([A-Za-z_][A-Za-z0-9_]*)[ \t]*\{", re.M)
SEED = 0
NOT_A_CLAIM = (
    "a stability record: what BioForge answers for each shipped design at a fixed seed, and which"
    " candidates tie with that answer. It is not a check that the answer is right; no field is an"
    " accuracy, a match to a publication or a probability, and BioForge states no probability"
    " (genomeos.organism.forge.design_certainty)"
)


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def design_programs(root: Path | None = None) -> list[str]:
    """Repository-relative paths of every `.bio` file under PROGRAM_ROOTS that holds a design block."""
    root = repo_root() if root is None else root
    found = []
    for top in PROGRAM_ROOTS:
        for path in sorted((root / top).rglob("*.bio")):
            if DESIGN_BLOCK.search(path.read_text(encoding="utf-8", errors="replace")):
                found.append(path.relative_to(root).as_posix())
    return sorted(found)


def results_for(program: str, seed: int = SEED) -> tuple[DesignResult, ...]:
    """Every design in one program, run once per process at `seed`; callers must not mutate the results."""
    return _results(str(program), int(seed))  # one cache key however the caller spells the arguments


@cache
def _results(program: str, seed: int) -> tuple[DesignResult, ...]:
    return tuple(run_designs(parse_file(str(repo_root() / program)), None, seed=seed))


def answer_record(program: str, result: DesignResult) -> dict[str, Any]:
    b = result.best
    answered = b is not None and b.feasible
    tied = (
        [
            c.label()
            for c in result.candidates
            if c.feasible and abs(c.loss - b.loss) <= 1e-12 and c.perturbations == b.perturbations
        ]
        if answered
        else []
    )
    return {
        "program": program,
        "design": result.design.name,
        "answer": b.label() if answered else None,
        "solved": result.solved,
        "loss": b.loss if answered else None,
        "tied_with_the_answer": [t for t in tied if t != b.label()],
        "tie_broken_by": "label order (forge.run_design's final sort)" if len(tied) > 1 else None,
        "evaluations": result.evaluations,
        "feasible_candidates": sum(1 for c in result.candidates if c.feasible),
        "varies_knobs": bool(result.design.vary),
    }


def answer_records(seed: int = SEED) -> list[dict[str, Any]]:
    """The stability record of every shipped design, in program then block order."""
    return [answer_record(p, r) for p in design_programs() for r in results_for(p, seed)]
