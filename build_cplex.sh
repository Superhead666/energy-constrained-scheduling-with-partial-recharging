#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

CPLEX_STUDIO_DIR="${CPLEX_STUDIO_DIR:-/Applications/CPLEX_Studio2211}"
CPLEX_ARCH="${CPLEX_ARCH:-arm64_osx}"
OUT_DIR="${SCRIPT_DIR}/bin"
mkdir -p "${OUT_DIR}"

g++ -std=c++17 -O3 -DNDEBUG -DWITH_CPLEX -DIL_STD \
  -I"${CPLEX_STUDIO_DIR}/cplex/include" \
  -I"${CPLEX_STUDIO_DIR}/concert/include" \
  "${SCRIPT_DIR}/main.cpp" \
  -L"${CPLEX_STUDIO_DIR}/cplex/lib/${CPLEX_ARCH}/static_pic" \
  -L"${CPLEX_STUDIO_DIR}/concert/lib/${CPLEX_ARCH}/static_pic" \
  -lilocplex -lconcert -lcplex -lm -lpthread -ldl \
  -o "${OUT_DIR}/caie_part1_final"

echo "${OUT_DIR}/caie_part1_final"
