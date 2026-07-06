#!/usr/bin/env bash
set -euo pipefail
AUX_PATH="${1:-testbench/ispd2005/adaptec1/adaptec1.aux}"
BUILD_DIR="${BUILD_DIR:-build}"
cmake -S . -B "$BUILD_DIR" ${LIMBO_ROOT:+-DLIMBO_ROOT="$LIMBO_ROOT"}
cmake --build "$BUILD_DIR" -j"$(nproc)"
"$BUILD_DIR/global_placer" --aux "$AUX_PATH"
