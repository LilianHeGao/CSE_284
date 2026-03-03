#!/usr/bin/env bash
set -euo pipefail

for cmd in bcftools plink2 gemma gcta64 python; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "Missing dependency: $cmd"
    exit 1
  fi
done

echo "All required commands found."
