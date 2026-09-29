#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Retired 2026-09-29. This script pushed dev, opened a pull request from dev into main, waited for CI
# and merged it. No session opens or merges a pull request, and main now moves through a branch
# frozen at a chosen sha that scripts/promote_main.sh --prepare pushes, in a pull request the owner
# opens and merges (CONTRIBUTING.md, "Release checks"). The file stays, refusing, because the
# checkout refuses file deletion; it runs nothing.
echo "promote.sh is retired: use scripts/promote_main.sh (CONTRIBUTING.md, \"Release checks\")" >&2
exit 1
