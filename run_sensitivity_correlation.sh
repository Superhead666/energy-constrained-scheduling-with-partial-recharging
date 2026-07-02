#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="${SCRIPT_DIR}/bin/caie_part1_final"

if [[ ! -x "${BIN}" ]]; then
  "${SCRIPT_DIR}/build_cplex.sh" >/dev/null
fi

N_EXACT="${N_EXACT:-50}"
N_ROBUST="${N_ROBUST:-100}"
TAU="${TAU:-0.5}"
R="${R:-0.4}"
GAMMA="${GAMMA:-4}"
DUE_REFERENCE_GAMMA="${DUE_REFERENCE_GAMMA:-4}"
DUE_REFERENCE_CORRELATION="${DUE_REFERENCE_CORRELATION:-independent}"
CORRELATIONS="${CORRELATIONS:-positive,independent,negative}"
REPS="${REPS:-20}"
TIME_EXACT="${TIME_EXACT:-600}"
TIME_HEURISTIC="${TIME_HEURISTIC:-30}"
WORKERS="${WORKERS:-4}"

echo "Running processing-energy correlation sensitivity"
echo "Exact grid: n=${N_EXACT}; tau=${TAU}; R=${R}; gamma=${GAMMA}; correlations=${CORRELATIONS}; reps=${REPS}; time=${TIME_EXACT}s"
echo "Robust grid: n=${N_ROBUST}; tau=${TAU}; R=${R}; gamma=${GAMMA}; correlations=${CORRELATIONS}; reps=${REPS}; time=${TIME_HEURISTIC}s"

N_EXACT="${N_EXACT}" \
N_ROBUST="${N_ROBUST}" \
TAU="${TAU}" \
R="${R}" \
GAMMA="${GAMMA}" \
DUE_REFERENCE_GAMMA="${DUE_REFERENCE_GAMMA}" \
DUE_REFERENCE_CORRELATION="${DUE_REFERENCE_CORRELATION}" \
CORRELATIONS="${CORRELATIONS}" \
REPS="${REPS}" \
TIME_EXACT="${TIME_EXACT}" \
TIME_HEURISTIC="${TIME_HEURISTIC}" \
WORKERS="${WORKERS}" \
"${SCRIPT_DIR}/run_sensitivity_correlation.py"
