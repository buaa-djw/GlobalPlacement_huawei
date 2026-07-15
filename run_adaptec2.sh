#!/usr/bin/env bash

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

AUX_FILE="${1:-${REPO_ROOT}/testbench/ispd2005/adaptec2/adaptec2.aux}"

EXECUTABLE="${REPO_ROOT}/build/placement_diagnostics"

if [[ ! -x "${EXECUTABLE}" ]]; then
    echo "Error: executable not found: ${EXECUTABLE}" >&2
    echo "Please build the project first." >&2
    exit 1
fi

if [[ ! -f "${AUX_FILE}" ]]; then
    echo "Error: benchmark not found: ${AUX_FILE}" >&2
    exit 1
fi

"${EXECUTABLE}" \
    --aux "${AUX_FILE}" \
    --bins 64 64 \
    --target-density 0.9 \
    --lambda 1.0 \
    --seed 20250715 \
    --samples 10