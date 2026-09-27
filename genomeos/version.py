# SPDX-License-Identifier: Apache-2.0
# Part of the BioLang engine (language, IR, VM, standard library); see LICENSING.md.
"""The version string, on the engine side of the licence split.

`bio` prints a version, and `bio` is Apache-2.0. If it read the version from
the package root it would import an AGPL module, which would make the
toolchain a work based on the application (LICENSING.md: the application may
import the engine, never the reverse). One line of engine-side code avoids
that, and it is also one less thing to untangle when BioLang is packaged on
its own, since the package root does not travel with it.

This string is the GenomeOS application's version. The packaged engine carries
its own (scripts/package_engine.py, ENGINE_VERSION), because the two are
released on different cadences.
"""

__version__ = "1.0.0"
