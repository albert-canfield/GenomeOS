# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The BioLang engine must not depend on the GenomeOS application.

This was an architectural rule (docs/ARCHITECTURE.md §10) before it was a
licensing one. Since 2026-09-11 the engine is Apache 2.0 and the application
is AGPL-3.0-or-later, so an engine file that imports an application module
would be a work based on AGPL code and could not honestly be distributed as
Apache (LICENSING.md, decision D40). That makes this test a licence check as
much as a design check, which is why it fails loudly rather than warning.

It also guards the roadmap's 2.0 milestone: the engine can only be packaged
separately for as long as these imports stay one-way.
"""

from __future__ import annotations

import ast
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

#: Everything the Apache-licensed half is made of.
ENGINE_PACKAGES = ("lang", "ir", "runtime", "std")
ENGINE_MODULES = ("coords", "version", "bio")

#: What an engine file is allowed to import from within the project.
ALLOWED = set(ENGINE_PACKAGES) | set(ENGINE_MODULES)


def engine_files() -> list[Path]:
    files = [p for pkg in ENGINE_PACKAGES for p in (ROOT / "genomeos" / pkg).rglob("*.py")]
    files += [ROOT / "genomeos" / f"{m}.py" for m in ENGINE_MODULES]
    return [p for p in files if p.exists() and "__pycache__" not in p.parts]


def imported_genomeos_names(path: Path) -> set[str]:
    """Every `genomeos.X` this file imports, at module level or inside a function."""
    names: set[str] = set()
    tree = ast.parse(path.read_text(), filename=str(path))
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module and node.module.startswith("genomeos"):
            parts = node.module.split(".")
            names.add(parts[1] if len(parts) > 1 else "genomeos")
        elif isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name.startswith("genomeos"):
                    parts = alias.name.split(".")
                    names.add(parts[1] if len(parts) > 1 else "genomeos")
    return names


def test_the_engine_is_a_real_set_of_files():
    files = engine_files()
    assert len(files) > 15, "the engine file list looks wrong; check ENGINE_PACKAGES"
    assert any(p.name == "parser.py" for p in files)
    assert any(p.name == "bio.py" for p in files)


def test_no_engine_file_imports_the_application():
    offenders: list[str] = []
    for path in engine_files():
        for name in imported_genomeos_names(path) - ALLOWED:
            rel = path.relative_to(ROOT)
            what = "the package root (genomeos/__init__.py)" if name == "genomeos" else f"genomeos.{name}"
            offenders.append(f"{rel} imports {what}")
    assert not offenders, (
        "the Apache-2.0 engine must not import the AGPL application "
        "(LICENSING.md, decision D40):\n  " + "\n  ".join(offenders)
    )


def test_the_package_root_is_not_an_engine_dependency():
    """`from genomeos import __version__` would pull in the AGPL package root."""
    for path in engine_files():
        assert "from genomeos import" not in path.read_text(), (
            f"{path.relative_to(ROOT)} reads the package root; import genomeos.version instead"
        )


def test_importing_the_toolchain_loads_no_application_module():
    """The strongest form of the rule: what actually gets loaded at run time."""
    code = (
        "import genomeos.bio, sys; "
        "print(' '.join(sorted(m for m in sys.modules if m.startswith('genomeos'))))"
    )
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, cwd=ROOT)
    assert out.returncode == 0, out.stderr
    loaded = out.stdout.split()
    assert loaded, "nothing was loaded; the subprocess probably failed"
    application = [m for m in loaded if m != "genomeos" and m.split(".")[1] not in ALLOWED]
    assert not application, f"`import genomeos.bio` pulled in application modules: {application}"
