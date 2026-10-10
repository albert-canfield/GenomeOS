# SPDX-License-Identifier: AGPL-3.0-or-later
"""Write data/results/forge_calibration_inventory.json: BioForge's calibration inventory (area H, 2026-10-10).

Applies the registration in genomeos/forge/calibration_inventory.py (committed in 5df76f2 before this
census was coded) to every design and experiment block under the CI `bio test` roots. Runs the four
shipped designs once (about a minute); poses no new design question to BioForge; fetches nothing.

    uv run python scripts/forge_calibration_inventory.py [--check]

`--check` writes nothing and exits 1 when the committed file differs from what the code computes now.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from genomeos.forge.calibration_inventory import build_inventory  # noqa: E402

OUT = ROOT / "data" / "results" / "forge_calibration_inventory.json"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true", help="compare with the committed file; write nothing")
    args = ap.parse_args()
    text = json.dumps(build_inventory(ROOT), indent=1, sort_keys=False) + "\n"
    if args.check:
        same = OUT.is_file() and OUT.read_text(encoding="utf-8") == text
        print("forge calibration inventory: " + ("matches the committed file" if same else "DIFFERS"))
        return 0 if same else 1
    OUT.write_text(text, encoding="utf-8")
    inv = json.loads(text)
    c = inv["counts"]
    print(
        f"wrote {OUT.relative_to(ROOT)}: {c['design_blocks']} design triples,"
        f" {c['experiment_blocks']} experiment blocks"
        f" ({c['expressible_as_inverse_design_questions']} expressible as design questions);"
        f" independent {c['independent_triples']} of a registered minimum {c['minimum_independent_triples']}."
        f" {inv['decision']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
