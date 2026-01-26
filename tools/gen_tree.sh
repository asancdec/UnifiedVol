#!/usr/bin/env bash
set -euo pipefail

ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT"

# If tree is not installed, skip silently
if ! command -v tree >/dev/null 2>&1; then
  echo "[tree] 'tree' not installed – skipping"
  exit 0
fi

tree -a -L 4 --dirsfirst \
  -I '.git|build|out|vcpkg|_deps|cmake-build-*|.cache|__pycache__|*.o|*.obj|*.so|*.a' \
  > /tmp/unifiedvol_tree.txt

{
  echo "# Repository tree"
  echo
  echo "_Auto-generated. Do not edit by hand._"
  echo
  echo '```text'
  cat /tmp/unifiedvol_tree.txt
  echo '```'
} > docs/TREE.md

rm -f /tmp/unifiedvol_tree.txt
