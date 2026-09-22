#!/usr/bin/env bash
# Install the pre-push check into this clone (see CONTRIBUTING.md).
set -euo pipefail
cd "$(dirname "$0")/.."
cp scripts/pre-push.sh .git/hooks/pre-push
chmod +x .git/hooks/pre-push
echo "pre-push hook installed"
