#!/usr/bin/env bash
set -euo pipefail

CFG=${1:-config/project.env}
source "$CFG"

if [[ -f "${CFG}.local" ]]; then
  source "${CFG}.local"
fi

check_tool() {
  local tool=$1
  local label=$2

  if [[ "$tool" == */* ]]; then
    if [[ ! -x "$tool" ]]; then
      echo "Missing dependency: $label ($tool)"
      exit 1
    fi
    return
  fi

  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "Missing dependency: $label"
    exit 1
  fi
}

check_tool "$PLINK2" plink2
check_tool python python

if [[ "$BCFTOOLS" == */* ]]; then
  if [[ ! -x "$BCFTOOLS" ]]; then
    echo "Optional dependency not found: bcftools ($BCFTOOLS)"
  fi
elif ! command -v "$BCFTOOLS" >/dev/null 2>&1; then
  echo "Optional dependency not found: bcftools"
fi

if [[ "$GEMMA" == */* ]]; then
  if [[ ! -x "$GEMMA" ]]; then
    echo "Optional dependency not found: gemma ($GEMMA)"
  fi
elif ! command -v "$GEMMA" >/dev/null 2>&1; then
  echo "Optional dependency not found: gemma"
fi

if [[ "$GCTA" == */* ]]; then
  if [[ ! -x "$GCTA" ]]; then
    echo "Optional dependency not found: gcta64 ($GCTA)"
  fi
elif ! command -v "$GCTA" >/dev/null 2>&1; then
  echo "Optional dependency not found: gcta64"
fi

echo "All required commands found."
